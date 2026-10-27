///|
/// System stability boundary and numerical stiffness analysis for RK4 integration.
///
/// Mathematical foundation:
/// For a linear time-invariant ODE dx/dt = A x with eigenvalues lambda_i = |lambda_i| * exp(j * theta_i),
/// Classical 4th-order Runge-Kutta (RK4) advances each eigenmode by the stability polynomial:
///   R(z) = 1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24,   where z = lambda_i * dt.
///
/// IMPORTANT: The commonly cited constant 2.78529 is ONLY the intercept of |R(z)| = 1 with the
/// negative real axis (theta = pi, damping ratio zeta = 1). For complex-conjugate poles in the
/// left-half plane (0 < zeta < 1), the boundary radius r_crit(theta) = |z_crit(theta)| dips below
/// 2.78529, reaching a minimum of r_min ≈ 2.56574 near zeta ≈ 0.75 (theta ≈ 138.6 deg).
/// Furthermore, if any pole has Re(lambda_i) > 0, the open-loop system is unstable for ALL dt > 0.
pub(all) struct SystemStabilityBound {
  max_pole_magnitude : Double
  min_pole_magnitude : Double
  stiffness_ratio : Double
  min_directional_radius : Double
  max_stable_dt : Double
  recommended_dt : Double
  is_stiff : Bool
  is_open_loop_stable : Bool
} derive(Eq)

pub extend SystemStabilityBound with Eq::{equal, not_equal}
pub extend SystemStabilityBound with Show::{output}

///|
/// Evaluate RK4 stability polynomial amplification factor |R(z)| for z = lambda * dt:
/// R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24.
/// A linear mode lambda is strictly contractive under RK4 if and only if |R(lambda * dt)| < 1.0.
pub fn rk4_amplification_factor(lambda : @num.Complex, dt : Double) -> Double {
  let z = lambda.mul_real(dt)
  let one = @num.Complex::new(1.0, 0.0)
  let z2 = (z * z).mul_real(1.0 / 2.0)
  let z3 = (z * z * z).mul_real(1.0 / 6.0)
  let z4 = (z * z * z * z).mul_real(1.0 / 24.0)
  let r = one + z + z2 + z3 + z4
  r.abs()
}

///|
/// Compute the exact directional stability boundary radius r_crit(theta) > 0 along the ray
/// z = r * (lambda / |lambda|) such that |R(r * lambda / |lambda|)| = 1.0.
/// Returns 0.0 if Re(lambda) > 0 (right-half plane pole is unstable for all r > 0).
pub fn rk4_directional_stability_radius(lambda : @num.Complex) -> Double {
  let mag = lambda.abs()
  if mag <= 1.0e-15 {
    return 0.0
  }
  // Strictly right-half plane pole: |R(r * u)| > 1 for small r > 0 -> unstable for all dt > 0
  if lambda.re > 1.0e-12 {
    return 0.0
  }
  // Purely imaginary pole: |R(j y)|^2 = 1 + y^6 (y^2 - 8) / 576 <= 1 iff 0 <= y <= 2*sqrt(2)
  if lambda.re.abs() <= 1.0e-12 {
    return 2.0 * (2.0 : Double).sqrt()
  }
  // Strictly left-half plane pole (Re(lambda) < 0):
  // Bisect on r in [0.1, 3.0] to find |R(r * u)| = 1.0 to machine precision.
  let unit = lambda.div_real(mag)
  let mut lo = 0.1
  let mut hi = 3.0
  for iter = 0; iter < 60; iter = iter + 1 {
    let mid = 0.5 * (lo + hi)
    let amp = rk4_amplification_factor(unit, mid)
    if amp <= 1.0 {
      lo = mid
    } else {
      hi = mid
    }
  }
  0.5 * (lo + hi)
}

///|
/// Compute the exact maximum stable step size dt_crit = r_crit(arg(lambda)) / |lambda|
/// for a single eigenvalue lambda under RK4 integration.
pub fn rk4_pole_max_stable_dt(lambda : @num.Complex) -> Double {
  let mag = lambda.abs()
  if mag <= 1.0e-15 {
    return 0.0
  }
  let r_crit = rk4_directional_stability_radius(lambda)
  r_crit / mag
}

///|
/// Verify whether a candidate step size dt is strictly inside the RK4 absolute stability region
/// for all poles of a transfer function.
pub fn verify_rk4_step_stability(
  tf : @tf.TransferFunction,
  dt : Double
) -> Bool {
  if dt <= 0.0 {
    return false
  }
  let poles = tf.poles()
  if poles.length() == 0 {
    return true
  }
  for pole in poles {
    if pole.re >= -1.0e-12 {
      return false
    }
    let amp = rk4_amplification_factor(pole, dt)
    if amp >= 1.0 {
      return false
    }
  }
  true
}

///|
/// Check exact 2D complex-plane stability boundaries, pole distribution, and stiffness ratio
/// for a given transfer function under Classical RK4 integration.
pub fn check_system_boundaries(tf : @tf.TransferFunction) -> SystemStabilityBound {
  let poles = tf.poles()
  let mut max_mag = 0.0
  let mut max_abs_re = 0.0
  let mut min_abs_re = 1.0e12
  let mut all_stable = true
  let mut min_dt_crit = 1.0e12
  let mut min_r_crit = 1.0e12

  for pole in poles {
    let re = pole.re
    // Strictly Hurwitz stable requires Re(pole) < 0 for all poles
    if re >= -1.0e-9 {
      all_stable = false
    }
    let mag = pole.abs()
    if mag > max_mag {
      max_mag = mag
    }
    let abs_re = re.abs()
    if abs_re > max_abs_re {
      max_abs_re = abs_re
    }
    if abs_re > 1.0e-9 && abs_re < min_abs_re {
      min_abs_re = abs_re
    }

    if re < -1.0e-9 && mag > 1.0e-12 {
      let r_crit = rk4_directional_stability_radius(pole)
      let dt_crit = r_crit / mag
      if r_crit < min_r_crit {
        min_r_crit = r_crit
      }
      if dt_crit < min_dt_crit {
        min_dt_crit = dt_crit
      }
    }
  }

  if min_abs_re > max_abs_re || min_abs_re >= 1.0e11 {
    min_abs_re = if max_abs_re > 1.0e-9 { max_abs_re } else { 1.0 }
  }
  let effective_max_re = if max_abs_re > 1.0e-9 { max_abs_re } else { 1.0 }
  let stiffness_ratio = effective_max_re / min_abs_re
  let is_stiff = stiffness_ratio > 100.0

  // If the system is open-loop unstable (Re(lambda) >= 0), no positive step size yields
  // an asymptotically stable open-loop trajectory: max_stable_dt is 0.0.
  let max_stable_dt = if all_stable && min_dt_crit < 1.0e11 {
    min_dt_crit
  } else {
    0.0
  }
  let min_directional_radius = if all_stable && min_r_crit < 1.0e11 {
    min_r_crit
  } else {
    0.0
  }

  // Recommended accuracy step size keeps dimensionless step |lambda_max| * dt <= 0.05
  let recommended_dt = if max_mag > 1.0e-9 {
    0.05 / max_mag
  } else {
    0.01
  }

  {
    max_pole_magnitude: max_mag,
    min_pole_magnitude: min_abs_re,
    stiffness_ratio,
    min_directional_radius,
    max_stable_dt,
    recommended_dt,
    is_stiff,
    is_open_loop_stable: all_stable,
  }
}

///|
/// Compute the spectral radius rho(Phi_CL(dt)) of the sampled-data closed-loop transition matrix
/// for a continuous plant controlled by a discrete PID controller at sampling/integration period dt.
///
/// Mathematical significance:
/// Even if the open-loop plant satisfies dt < r_crit / |lambda_plant| AND the continuous-time
/// closed-loop transfer function is Hurwitz stable, the sampled-data ZOH + discrete PID loop
/// is asymptotically stable if and only if rho(Phi_CL(dt)) < 1.0.
/// Uses Gelfand's formula rho(M) = lim_{k->inf} ||M^(2^k)||_inf^(1/2^k) with 7 squarings (M^128).
pub fn closed_loop_pid_spectral_radius(
  plant : @tf.TransferFunction,
  kp : Double,
  ki : Double,
  kd : Double,
  dt : Double
) -> Result[Double, String] {
  if dt <= 0.0 {
    return Err("Step size dt must be strictly positive")
  }
  let model = match StateSpaceModel::from_tf(plant) {
    Ok(m) => m
    Err(e) => return Err(e)
  }
  let n = model.order
  // Augmented state z = [x_0, ..., x_{n-1}, integral, prev_error] of dimension m = n + 2
  let m = n + 2
  let phi = @num.Matrix::zeros(m, m)

  // Construct exact transition matrix column by column using standard basis vectors e_col
  for col = 0; col < m; col = col + 1 {
    let x = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      if i == col {
        x[i] = 1.0
      }
    }
    let integral_prev = if col == n { 1.0 } else { 0.0 }
    let err_prev = if col == n + 1 { 1.0 } else { 0.0 }

    let y = model.compute_output(x, 0.0)
    let err = -y
    let integral_next = integral_prev + err * dt
    let deriv = (err - err_prev) / dt
    let u = kp * err + ki * integral_next + kd * deriv
    let x_next = rk4_step(model, x, u, dt)

    for row = 0; row < n; row = row + 1 {
      phi.set(row, col, x_next[row])
    }
    phi.set(n, col, integral_next)
    phi.set(n + 1, col, err)
  }

  // Normalized Gelfand repeated squaring over 14 steps (Phi^(2^14) = Phi^16384):
  // Avoids overflow/underflow while eliminating initial non-normal transient prefactors.
  let mut mat = phi
  let mut log_rho = 0.0
  let mut weight = 0.5
  for step = 0; step < 14; step = step + 1 {
    mat = mat * mat
    let mut norm_inf = 0.0
    for r = 0; r < m; r = r + 1 {
      let mut row_sum = 0.0
      for c = 0; c < m; c = c + 1 {
        row_sum = row_sum + mat.get(r, c).abs()
      }
      if row_sum > norm_inf {
        norm_inf = row_sum
      }
    }
    if norm_inf <= 1.0e-300 {
      return Ok(0.0)
    }
    log_rho = log_rho + weight * @math.ln(norm_inf)
    mat = mat.scale(1.0 / norm_inf)
    weight = weight * 0.5
  }
  Ok(@math.exp(log_rho))
}

///|
pub fn SystemStabilityBound::to_string(self : SystemStabilityBound) -> String {
  "SystemStabilityBound(max_pole: \{self.max_pole_magnitude}, stiffness_ratio: \{self.stiffness_ratio}, min_r_crit: \{self.min_directional_radius}, max_dt: \{self.max_stable_dt}, recommended_dt: \{self.recommended_dt}, stiff: \{self.is_stiff}, stable: \{self.is_open_loop_stable})"
}

///|
pub impl Show for SystemStabilityBound with output(self, logger) {
  logger.write_string(self.to_string())
}
