///|
/// Result of a time-domain simulation containing sampled time points, outputs, and inputs.
pub(all) struct SimulationResult {
  times : Array[Double]
  outputs : Array[Double]
  inputs : Array[Double]
} derive(Eq)

pub extend SimulationResult with Eq::{equal, not_equal}
pub extend SimulationResult with Show::{output}

///|
/// Internal state space model in controllable canonical form:
/// dx/dt = A*x + B*u
/// y = C*x + D*u
priv struct StateSpaceModel {
  order : Int
  a_coeffs : Array[Double] // [a_0, a_1, ..., a_{n-1}]
  b_coeffs : Array[Double] // [b_0, b_1, ..., b_{n-1}]
  d_gain : Double
}

///|
/// Convert transfer function to controllable canonical state space realization.
fn StateSpaceModel::from_tf(
  tf : @tf.TransferFunction
) -> Result[StateSpaceModel, String] {
  let den_deg = tf.den.degree()
  let num_deg = tf.num.degree()

  if den_deg == 0 {
    return Err("Cannot simulate system of order 0 with differential equations")
  }
  if num_deg > den_deg {
    return Err("Transfer function is improper: order(num) > order(den)")
  }

  let n = den_deg
  let a_coeffs = Array::make(n, 0.0)
  let b_coeffs = Array::make(n, 0.0)

  // Denominator is monic [1.0, a_{n-1}, a_{n-2}, ..., a_0]
  // In coeffs array: den.coeffs[0] = 1.0, den.coeffs[1] = a_{n-1}, ..., den.coeffs[n] = a_0
  for i = 0; i < n; i = i + 1 {
    let power = i
    let idx_in_den = n - power
    a_coeffs[i] = tf.den.coeffs[idx_in_den]
  }

  // Handle direct feedthrough d if num_deg == den_deg
  let mut d_gain = 0.0
  let num_coeffs_len = tf.num.coeffs.length()
  if num_deg == den_deg {
    d_gain = tf.num.coeffs[0]
    // Remainder numerator = num - d * den
    for i = 0; i < n; i = i + 1 {
      let power = i
      let idx_in_num = num_coeffs_len - 1 - power
      let num_c = if idx_in_num >= 0 { tf.num.coeffs[idx_in_num] } else { 0.0 }
      let den_c = a_coeffs[i]
      b_coeffs[i] = num_c - d_gain * den_c
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      let power = i
      let idx_in_num = num_coeffs_len - 1 - power
      let num_c = if idx_in_num >= 0 { tf.num.coeffs[idx_in_num] } else { 0.0 }
      b_coeffs[i] = num_c
    }
  }

  Ok({ order: n, a_coeffs, b_coeffs, d_gain })
}

///|
/// State derivative: dx/dt = f(x, u)
fn StateSpaceModel::derivative(
  self : StateSpaceModel,
  x : Array[Double],
  u : Double
) -> Array[Double] {
  let n = self.order
  let dx = Array::make(n, 0.0)
  for i = 0; i < n - 1; i = i + 1 {
    dx[i] = x[i + 1]
  }
  let mut sum_ax = 0.0
  for i = 0; i < n; i = i + 1 {
    sum_ax = sum_ax + self.a_coeffs[i] * x[i]
  }
  dx[n - 1] = -sum_ax + u
  dx
}

///|
/// System output: y = C*x + D*u
fn StateSpaceModel::compute_output(
  self : StateSpaceModel,
  x : Array[Double],
  u : Double
) -> Double {
  let mut y = self.d_gain * u
  for i = 0; i < self.order; i = i + 1 {
    y = y + self.b_coeffs[i] * x[i]
  }
  y
}

///|
/// Single step of Runge-Kutta 4 (RK4) integration:
/// advances state vector x by time dt under control input u.
fn rk4_step(
  model : StateSpaceModel,
  x : Array[Double],
  u : Double,
  dt : Double
) -> Array[Double] {
  let n = model.order

  // k1 = f(x, u)
  let k1 = model.derivative(x, u)

  // x_k1 = x + 0.5 * dt * k1
  let x_k1 = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    x_k1[i] = x[i] + 0.5 * dt * k1[i]
  }
  let k2 = model.derivative(x_k1, u)

  // x_k2 = x + 0.5 * dt * k2
  let x_k2 = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    x_k2[i] = x[i] + 0.5 * dt * k2[i]
  }
  let k3 = model.derivative(x_k2, u)

  // x_k3 = x + dt * k3
  let x_k3 = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    x_k3[i] = x[i] + dt * k3[i]
  }
  let k4 = model.derivative(x_k3, u)

  // x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
  let x_next = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    x_next[i] = x[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
  }
  x_next
}

///|
/// Simulate open-loop unit step response of any proper transfer function using RK4.
pub fn simulate_step(
  tf : @tf.TransferFunction,
  duration : Double,
  dt~ : Double = 0.001
) -> Result[SimulationResult, String] {
  if duration <= 0.0 {
    return Err("Duration must be strictly positive")
  }
  if dt <= 0.0 {
    return Err("Time step dt must be strictly positive")
  }

  let model = match StateSpaceModel::from_tf(tf) {
    Ok(m) => m
    Err(e) => return Err(e)
  }
  let num_steps = (duration / dt).to_int() + 1

  let times = Array(capacity=num_steps)
  let outputs = Array(capacity=num_steps)
  let inputs = Array(capacity=num_steps)

  let mut state = Array::make(model.order, 0.0)
  let u = 1.0

  for step = 0; step < num_steps; step = step + 1 {
    let t = step.to_double() * dt
    let y = model.compute_output(state, u)
    times.push(t)
    outputs.push(y)
    inputs.push(u)
    state = rk4_step(model, state, u, dt)
  }

  Ok({ times, outputs, inputs })
}

///|
/// Simulate closed-loop feedback control with a plant transfer function and a discrete PID controller.
pub fn simulate_closed_loop_pid(
  plant : @tf.TransferFunction,
  pid : @pid.PIDController,
  duration : Double,
  setpoint~ : Double = 1.0,
  dt~ : Double = 0.001
) -> Result[SimulationResult, String] {
  if duration <= 0.0 {
    return Err("Duration must be strictly positive")
  }
  if dt <= 0.0 {
    return Err("Time step dt must be strictly positive")
  }

  let model = match StateSpaceModel::from_tf(plant) {
    Ok(m) => m
    Err(e) => return Err(e)
  }
  let num_steps = (duration / dt).to_int() + 1

  let times = Array(capacity=num_steps)
  let outputs = Array(capacity=num_steps)
  let inputs = Array(capacity=num_steps)

  let mut state = Array::make(model.order, 0.0)
  pid.reset()

  for step = 0; step < num_steps; step = step + 1 {
    let t = step.to_double() * dt
    let y = model.compute_output(state, 0.0) // current output before new control
    let u = pid.update(setpoint, y, dt)

    times.push(t)
    outputs.push(y)
    inputs.push(u)

    state = rk4_step(model, state, u, dt)
  }

  Ok({ times, outputs, inputs })
}

///|
/// Compute maximum absolute error between simulation outputs and 1st-order analytical response.
pub fn error_against_first_order_analytic(
  k : Double,
  t_const : Double,
  sim_result : SimulationResult
) -> Result[Double, String] {
  let mut max_err = 0.0
  let len = sim_result.times.length()
  for i = 0; i < len; i = i + 1 {
    let t = sim_result.times[i]
    let y_sim = sim_result.outputs[i]
    let y_exact = match @analysis.step_response_first_order(k, t_const, t) {
      Ok(val) => val
      Err(err) => return Err(err)
    }
    let err = (y_sim - y_exact).abs()
    if err > max_err {
      max_err = err
    }
  }
  Ok(max_err)
}

///|
/// Compute maximum absolute error between simulation outputs and 2nd-order analytical response.
pub fn error_against_second_order_analytic(
  omega_n : Double,
  zeta : Double,
  sim_result : SimulationResult
) -> Result[Double, String] {
  let mut max_err = 0.0
  let len = sim_result.times.length()
  for i = 0; i < len; i = i + 1 {
    let t = sim_result.times[i]
    let y_sim = sim_result.outputs[i]
    let y_exact = match @analysis.step_response_second_order(omega_n, zeta, t) {
      Ok(val) => val
      Err(err) => return Err(err)
    }
    let err = (y_sim - y_exact).abs()
    if err > max_err {
      max_err = err
    }
  }
  Ok(max_err)
}

///|
/// Quantitative time-domain performance metrics and integral error indices measured directly
/// from a simulated trajectory (supports both linear and nonlinear saturated/delayed loops).
pub(all) struct TrajectoryMetrics {
  overshoot_percent : Double
  peak_value : Double
  peak_time : Double
  rise_time_10_90 : Double
  settling_time_2pct : Double
  steady_state_error : Double
  iae : Double
  ise : Double
  itae : Double
  control_effort_l2 : Double
} derive(Eq)

pub extend TrajectoryMetrics with Eq::{equal, not_equal}
pub extend TrajectoryMetrics with Show::{output}

///|
/// Measure engineering performance metrics and integral indices (IAE, ISE, ITAE) directly
/// from any simulated time-domain trajectory.
pub fn measure_trajectory_metrics(
  sim : SimulationResult,
  setpoint~ : Double = 1.0
) -> TrajectoryMetrics {
  let n = sim.times.length()
  if n == 0 {
    return {
      overshoot_percent: 0.0,
      peak_value: 0.0,
      peak_time: 0.0,
      rise_time_10_90: 0.0,
      settling_time_2pct: 0.0,
      steady_state_error: 0.0,
      iae: 0.0,
      ise: 0.0,
      itae: 0.0,
      control_effort_l2: 0.0,
    }
  }

  let mut peak_val = sim.outputs[0]
  let mut peak_t = sim.times[0]
  let mut t10 = -1.0
  let mut t90 = -1.0
  let target_10 = 0.1 * setpoint
  let target_90 = 0.9 * setpoint
  let band = 0.02 * (if setpoint.abs() > 1.0e-9 { setpoint.abs() } else { 1.0 })
  let mut last_outside_t = 0.0

  let mut iae = 0.0
  let mut ise = 0.0
  let mut itae = 0.0
  let mut u_l2 = 0.0

  for i = 0; i < n; i = i + 1 {
    let t = sim.times[i]
    let y = sim.outputs[i]
    let u = sim.inputs[i]
    let err = setpoint - y
    let abs_err = err.abs()

    if y > peak_val {
      peak_val = y
      peak_t = t
    }
    if t10 < 0.0 && y >= target_10 {
      t10 = t
    }
    if t90 < 0.0 && y >= target_90 {
      t90 = t
    }
    if abs_err > band {
      last_outside_t = t
    }

    if i > 0 {
      let dt = t - sim.times[i - 1]
      iae = iae + abs_err * dt
      ise = ise + err * err * dt
      itae = itae + t * abs_err * dt
      u_l2 = u_l2 + u * u * dt
    }
  }

  let overshoot_pct = if setpoint.abs() > 1.0e-9 && peak_val > setpoint {
    (peak_val - setpoint) / setpoint.abs() * 100.0
  } else {
    0.0
  }
  let rise_t = if t10 >= 0.0 && t90 >= t10 { t90 - t10 } else { 0.0 }
  let final_err = (setpoint - sim.outputs[n - 1]).abs()

  {
    overshoot_percent: overshoot_pct,
    peak_value: peak_val,
    peak_time: peak_t,
    rise_time_10_90: rise_t,
    settling_time_2pct: last_outside_t,
    steady_state_error: final_err,
    iae,
    ise,
    itae,
    control_effort_l2: u_l2,
  }
}

///|
/// Simulate a hybrid nonlinear closed-loop control system with:
/// - Continuous-time plant ODE integrated via multi-rate substeps (Zero-Order Hold)
/// - Actuator saturation [u_min, u_max] and slew-rate limit |du/dt| <= max_slew_rate
/// - Back-calculation anti-windup (tracking gain kb) and clamping
/// - Pure transport dead-time delay exp(-transport_delay * s) via exact FIFO ring buffer
/// - External load disturbance injection d(t) at t >= disturbance_time
///
/// Note: None of these nonlinear/hybrid/delayed phenomena can be simulated by symbolic
/// Laplace partial-fraction expansion packages such as CAIMEOX/symbit/symphysics/control.
pub fn simulate_closed_loop_nonlinear(
  plant : @tf.TransferFunction,
  kp : Double,
  ki : Double,
  kd : Double,
  duration : Double,
  setpoint~ : Double = 1.0,
  dt~ : Double = 0.005,
  substeps~ : Int = 4,
  u_min? : Double,
  u_max? : Double,
  max_slew_rate? : Double,
  anti_windup_kb~ : Double = 0.0,
  transport_delay~ : Double = 0.0,
  disturbance_time~ : Double = -1.0,
  disturbance_magnitude~ : Double = 0.0
) -> Result[SimulationResult, String] {
  if duration <= 0.0 || dt <= 0.0 {
    return Err("Duration and time step dt must be strictly positive")
  }
  let model = match StateSpaceModel::from_tf(plant) {
    Ok(m) => m
    Err(e) => return Err(e)
  }

  let m_sub = if substeps > 1 { substeps } else { 1 }
  let dt_sub = dt / m_sub.to_double()
  let num_steps = (duration / dt).to_int() + 1

  // Delay buffer in controller sample steps
  let delay_steps = if transport_delay > 0.0 {
    (transport_delay / dt + 0.5).to_int()
  } else {
    0
  }
  let delay_buf = Array::make(delay_steps + 1, 0.0)
  let mut delay_idx = 0

  let times = Array(capacity=num_steps)
  let outputs = Array(capacity=num_steps)
  let inputs = Array(capacity=num_steps)

  let mut state = Array::make(model.order, 0.0)
  let mut integral = 0.0
  let mut prev_err = 0.0
  let mut prev_u = 0.0
  let mut initialized = false

  for step = 0; step < num_steps; step = step + 1 {
    let t = step.to_double() * dt
    let y = model.compute_output(state, 0.0)
    let err = setpoint - y

    let deriv = if initialized {
      (err - prev_err) / dt
    } else {
      initialized = true
      0.0
    }
    prev_err = err

    // Unconstrained PID output before saturation
    let u_raw = kp * err + ki * integral + kd * deriv
    let mut u_cmd = u_raw

    // Apply actuator amplitude saturation
    match u_min {
      Some(lo) => if u_cmd < lo { u_cmd = lo }
      None => ()
    }
    match u_max {
      Some(hi) => if u_cmd > hi { u_cmd = hi }
      None => ()
    }

    // Apply actuator slew-rate limit |u_cmd - prev_u| <= max_slew_rate * dt
    match max_slew_rate {
      Some(rate) => if step > 0 {
        let max_du = rate.abs() * dt
        if u_cmd - prev_u > max_du {
          u_cmd = prev_u + max_du
        } else if prev_u - u_cmd > max_du {
          u_cmd = prev_u - max_du
        }
      }
      None => ()
    }
    prev_u = u_cmd

    // Update integrator with back-calculation anti-windup:
    // dI/dt = err + (kb / ki) * (u_cmd - u_raw) when ki != 0
    if ki.abs() > 1.0e-12 {
      let aw_correction = anti_windup_kb * (u_cmd - u_raw) / ki
      integral = integral + (err + aw_correction) * dt
    }

    // Transport delay FIFO queue
    let u_delayed = if delay_steps > 0 {
      let out_u = delay_buf[delay_idx]
      delay_buf[delay_idx] = u_cmd
      delay_idx = (delay_idx + 1) % delay_steps
      out_u
    } else {
      u_cmd
    }

    // External load disturbance injection
    let d_val = if disturbance_time >= 0.0 && t >= disturbance_time {
      disturbance_magnitude
    } else {
      0.0
    }
    let u_plant = u_delayed + d_val

    times.push(t)
    outputs.push(y)
    inputs.push(u_cmd)

    // Multi-rate continuous plant integration under Zero-Order Hold (ZOH)
    for sub = 0; sub < m_sub; sub = sub + 1 {
      state = rk4_step(model, state, u_plant, dt_sub)
    }
  }

  Ok({ times, outputs, inputs })
}

///|
pub fn TrajectoryMetrics::to_string(self : TrajectoryMetrics) -> String {
  "TrajectoryMetrics(Mp: \{self.overshoot_percent}%, tp: \{self.peak_time}s, tr: \{self.rise_time_10_90}s, ts: \{self.settling_time_2pct}s, e_ss: \{self.steady_state_error}, IAE: \{self.iae}, ISE: \{self.ise}, ITAE: \{self.itae})"
}

///|
pub impl Show for TrajectoryMetrics with output(self, logger) {
  logger.write_string(self.to_string())
}

///|
pub fn SimulationResult::to_string(self : SimulationResult) -> String {
  let len = self.times.length()
  if len == 0 {
    return "SimulationResult(empty)"
  }
  let final_t = self.times[len - 1]
  let final_y = self.outputs[len - 1]
  "SimulationResult(points: \{len}, t_final: \{final_t}s, y_final: \{final_y})"
}

///|
pub impl Show for SimulationResult with output(self, logger) {
  logger.write_string(self.to_string())
}
