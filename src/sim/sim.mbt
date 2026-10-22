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
