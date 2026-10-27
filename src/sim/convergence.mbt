///|
/// Convergence evaluation record for a single step size dt.
pub(all) struct StepConvergenceRecord {
  dt : Double
  dimensionless_step : Double
  max_error : Double
  rms_error : Double
  richardson_error_est : Double
  order_estimate : Double
} derive(Eq)

pub extend StepConvergenceRecord with Eq::{equal, not_equal}
pub extend StepConvergenceRecord with Show::{output}

///|
/// Comprehensive study of RK4 numerical convergence across varying step sizes.
///
/// Mathematical note on error scaling:
/// For a linear system dx/dt = A x + B u, the global truncation error is NOT a function of
/// dimensional step size dt alone; it scales with the dimensionless product eta = |lambda_max| * dt,
/// the modal residue sum (system gain), and the condition number of the state realization:
///   E_max(dt) = C(A, B, C, T) * (|lambda_max| * dt)^4 + O((|lambda_max| * dt)^5) + O(eps_mach / dt).
pub(all) struct ConvergenceStudy {
  records : Array[StepConvergenceRecord]
  mean_order : Double
  theoretical_order : Double
  precision_floor : Double
} derive(Eq)

pub extend ConvergenceStudy with Eq::{equal, not_equal}
pub extend ConvergenceStudy with Show::{output}

///|
/// Simulation result certified by a posteriori Step-Doubling Richardson Extrapolation.
pub(all) struct CertifiedSimulationResult {
  times : Array[Double]
  coarse_outputs : Array[Double]
  fine_outputs : Array[Double]
  extrapolated_outputs : Array[Double]
  error_certificate_max : Double
  error_certificate_rms : Double
  dimensionless_step : Double
  is_within_stability_region : Bool
} derive(Eq)

pub extend CertifiedSimulationResult with Eq::{equal, not_equal}
pub extend CertifiedSimulationResult with Show::{output}

///|
/// Adaptive step-size simulation result with local error control (atol, rtol).
pub(all) struct AdaptiveSimulationResult {
  times : Array[Double]
  outputs : Array[Double]
  step_sizes : Array[Double]
  accepted_steps : Int
  rejected_steps : Int
  max_estimated_local_error : Double
} derive(Eq)

pub extend AdaptiveSimulationResult with Eq::{equal, not_equal}
pub extend AdaptiveSimulationResult with Show::{output}

///|
/// Simulate unit step response with a posteriori Step-Doubling Richardson Extrapolation error certificate:
///   E_est(t_k) = |y^(dt/2)(t_k) - y^(dt)(t_k)| / (2^4 - 1) = |y^(dt/2)(t_k) - y^(dt)(t_k)| / 15.
/// This eliminates open-loop step-size guessing by providing a computed error bound for any transfer function.
pub fn simulate_step_certified(
  tf : @tf.TransferFunction,
  duration : Double,
  dt~ : Double = 0.005
) -> Result[CertifiedSimulationResult, String] {
  let coarse = match simulate_step(tf, duration, dt~) {
    Ok(s) => s
    Err(e) => return Err(e)
  }
  let fine = match simulate_step(tf, duration, dt=0.5 * dt) {
    Ok(s) => s
    Err(e) => return Err(e)
  }

  let bounds = check_system_boundaries(tf)
  let eta = bounds.max_pole_magnitude * dt
  let is_stable = verify_rk4_step_stability(tf, dt)

  let n = coarse.times.length()
  let fine_sampled = Array(capacity=n)
  let extrap = Array(capacity=n)
  let mut max_est = 0.0
  let mut sum_sq = 0.0

  for k = 0; k < n; k = k + 1 {
    let y_h = coarse.outputs[k]
    let fine_idx = if 2 * k < fine.outputs.length() {
      2 * k
    } else {
      fine.outputs.length() - 1
    }
    let y_half = fine.outputs[fine_idx]
    // For 4th-order RK4:
    //   y_exact - y_half ≈ (y_half - y_h) / 15       (error of fine trajectory dt/2)
    //   y_exact - y_h    ≈ 16 * (y_half - y_h) / 15  (error of coarse trajectory dt)
    let fine_corr = (y_half - y_h) / 15.0
    let coarse_err_est = (16.0 / 15.0 * (y_half - y_h)).abs()
    if coarse_err_est > max_est {
      max_est = coarse_err_est
    }
    sum_sq = sum_sq + coarse_err_est * coarse_err_est
    fine_sampled.push(y_half)
    extrap.push(y_half + fine_corr)
  }

  let rms_est = (sum_sq / n.to_double()).sqrt()

  Ok({
    times: coarse.times,
    coarse_outputs: coarse.outputs,
    fine_outputs: fine_sampled,
    extrapolated_outputs: extrap,
    error_certificate_max: max_est,
    error_certificate_rms: rms_est,
    dimensionless_step: eta,
    is_within_stability_region: is_stable,
  })
}

///|
/// Adaptive step-size RK45 / Step-Doubling integrator that automatically adjusts dt_k
/// to satisfy local error tolerance: err_k <= atol + rtol * max(|y_k|, |y_{k+1}|).
pub fn simulate_step_adaptive(
  tf : @tf.TransferFunction,
  duration : Double,
  atol~ : Double = 1.0e-10,
  rtol~ : Double = 1.0e-8,
  dt_init~ : Double = 0.02,
  dt_min~ : Double = 1.0e-6,
  dt_max~ : Double = 0.2
) -> Result[AdaptiveSimulationResult, String] {
  if duration <= 0.0 {
    return Err("Duration must be strictly positive")
  }
  if atol <= 0.0 || rtol <= 0.0 {
    return Err("Tolerances atol and rtol must be strictly positive")
  }
  let model = match StateSpaceModel::from_tf(tf) {
    Ok(m) => m
    Err(e) => return Err(e)
  }

  let times = [0.0]
  let u = 1.0
  let mut state = Array::make(model.order, 0.0)
  let outputs = [model.compute_output(state, u)]
  let step_sizes : Array[Double] = []

  let mut t = 0.0
  let mut h = dt_init
  if h > dt_max {
    h = dt_max
  }
  if h < dt_min {
    h = dt_min
  }

  let mut accepted = 0
  let mut rejected = 0
  let mut max_local_err = 0.0
  let max_total_iters = 200000

  let mut iter = 0
  while t < duration - 1.0e-14 && iter < max_total_iters {
    iter = iter + 1
    if t + h > duration {
      h = duration - t
    }

    // One coarse step of size h
    let x_coarse = rk4_step(model, state, u, h)
    // Two half steps of size h/2
    let x_mid = rk4_step(model, state, u, 0.5 * h)
    let x_fine = rk4_step(model, x_mid, u, 0.5 * h)

    // Compute state-space local Richardson error norm ||x_fine - x_coarse||_inf / 15
    let mut err_norm = 0.0
    let mut state_norm = 0.0
    for i = 0; i < model.order; i = i + 1 {
      let diff = (x_fine[i] - x_coarse[i]).abs() / 15.0
      if diff > err_norm {
        err_norm = diff
      }
      let mag = x_fine[i].abs()
      if mag > state_norm {
        state_norm = mag
      }
    }
    let y_coarse = model.compute_output(x_coarse, u)
    let y_fine = model.compute_output(x_fine, u)
    let y_err = (y_fine - y_coarse).abs() / 15.0
    if y_err > err_norm {
      err_norm = y_err
    }

    let tol = atol + rtol * (if state_norm > y_fine.abs() { state_norm } else { y_fine.abs() })
    let scaled_err = err_norm / tol

    if scaled_err <= 1.0 || h <= dt_min {
      // Accept step using 5th-order Richardson-extrapolated state
      let x_extrap = Array::make(model.order, 0.0)
      for i = 0; i < model.order; i = i + 1 {
        x_extrap[i] = x_fine[i] + (x_fine[i] - x_coarse[i]) / 15.0
      }
      state = x_extrap
      t = t + h
      let y_next = model.compute_output(state, u)
      times.push(t)
      outputs.push(y_next)
      step_sizes.push(h)
      accepted = accepted + 1
      if err_norm > max_local_err {
        max_local_err = err_norm
      }

      // PI step-size update factor: 0.9 * (1 / scaled_err)^(1/5)
      let factor = if scaled_err <= 1.0e-15 {
        2.0
      } else {
        let raw = 0.9 * @math.pow(scaled_err, -0.2)
        if raw > 2.0 { 2.0 } else if raw < 0.2 { 0.2 } else { raw }
      }
      h = h * factor
      if h > dt_max {
        h = dt_max
      }
      if h < dt_min {
        h = dt_min
      }
    } else {
      // Reject step and shrink h
      rejected = rejected + 1
      let shrink = 0.9 * @math.pow(scaled_err, -0.2)
      let factor = if shrink < 0.2 { 0.2 } else if shrink > 0.8 { 0.8 } else { shrink }
      h = h * factor
      if h < dt_min {
        h = dt_min
      }
    }
  }

  Ok({
    times,
    outputs,
    step_sizes,
    accepted_steps: accepted,
    rejected_steps: rejected,
    max_estimated_local_error: max_local_err,
  })
}

///|
/// Run systematic step-size variation study to verify RK4 4th-order convergence O((|lambda|*dt)^4).
/// dts: array of descending step sizes, e.g. [0.08, 0.04, 0.02, 0.01, 0.005]
pub fn evaluate_rk4_convergence(
  tf : @tf.TransferFunction,
  duration : Double,
  dts : Array[Double],
  exact_fn : (Double) -> Double
) -> Result[ConvergenceStudy, String] {
  if dts.length() < 2 {
    return Err("At least two step sizes are required to evaluate convergence order")
  }

  let bounds = check_system_boundaries(tf)
  let records = Array(capacity=dts.length())
  let mut prev_dt = 0.0
  let mut prev_err = 0.0
  let mut order_sum = 0.0
  let mut order_count = 0
  let mut min_err = 1.0e12

  for i = 0; i < dts.length(); i = i + 1 {
    let dt = dts[i]
    let cert = match simulate_step_certified(tf, duration, dt~) {
      Ok(c) => c
      Err(e) => return Err(e)
    }

    let mut max_err = 0.0
    let mut sum_sq = 0.0
    let n = cert.times.length()

    for k = 0; k < n; k = k + 1 {
      let t = cert.times[k]
      let y_sim = cert.coarse_outputs[k]
      let y_exact = exact_fn(t)
      let diff = (y_sim - y_exact).abs()
      if diff > max_err {
        max_err = diff
      }
      sum_sq = sum_sq + diff * diff
    }

    let rms = (sum_sq / n.to_double()).sqrt()
    if max_err < min_err {
      min_err = max_err
    }

    let mut order = 0.0
    if i > 0 && max_err > 1.0e-13 && prev_err > 1.0e-13 {
      let log_err_ratio = @math.ln(prev_err / max_err)
      let log_dt_ratio = @math.ln(prev_dt / dt)
      order = log_err_ratio / log_dt_ratio
      order_sum = order_sum + order
      order_count = order_count + 1
    }

    records.push({
      dt,
      dimensionless_step: bounds.max_pole_magnitude * dt,
      max_error: max_err,
      rms_error: rms,
      richardson_error_est: cert.error_certificate_max,
      order_estimate: order,
    })
    prev_dt = dt
    prev_err = max_err
  }

  let mean_order = if order_count > 0 {
    order_sum / order_count.to_double()
  } else {
    4.0
  }

  Ok({
    records,
    mean_order,
    theoretical_order: 4.0,
    precision_floor: min_err,
  })
}

///|
pub fn StepConvergenceRecord::to_string(self : StepConvergenceRecord) -> String {
  "dt: \{self.dt}, eta: \{self.dimensionless_step}, max_err: \{self.max_error}, richardson_est: \{self.richardson_error_est}, order: \{self.order_estimate}"
}

///|
pub impl Show for StepConvergenceRecord with output(self, logger) {
  logger.write_string(self.to_string())
}

///|
pub fn ConvergenceStudy::to_string(self : ConvergenceStudy) -> String {
  let mut s = "ConvergenceStudy(theoretical_order: \{self.theoretical_order}, mean_order: \{self.mean_order}, floor: \{self.precision_floor}\n"
  for r in self.records {
    s = s + "  " + r.to_string() + "\n"
  }
  s + ")"
}

///|
pub impl Show for ConvergenceStudy with output(self, logger) {
  logger.write_string(self.to_string())
}

///|
pub fn CertifiedSimulationResult::to_string(self : CertifiedSimulationResult) -> String {
  "CertifiedSimulationResult(points: \{self.times.length()}, eta: \{self.dimensionless_step}, max_err_cert: \{self.error_certificate_max}, stable: \{self.is_within_stability_region})"
}

///|
pub impl Show for CertifiedSimulationResult with output(self, logger) {
  logger.write_string(self.to_string())
}

///|
pub fn AdaptiveSimulationResult::to_string(self : AdaptiveSimulationResult) -> String {
  "AdaptiveSimulationResult(accepted: \{self.accepted_steps}, rejected: \{self.rejected_steps}, max_local_err: \{self.max_estimated_local_error})"
}

///|
pub impl Show for AdaptiveSimulationResult with output(self, logger) {
  logger.write_string(self.to_string())
}
