test "RK4 first-order step response matches analytical solution to < 1e-6" {
  // G(s) = 2.0 / (0.5s + 1.0) -> K = 2.0, T = 0.5
  let tf = @tf.TransferFunction::from_coeffs([2.0], [0.5, 1.0]).unwrap()
  let dt = 0.0005
  let sim = simulate_step(tf, 2.0, dt~).unwrap()

  // Check at t = 0.5, 1.0, 1.5, 2.0
  let test_times = [0.5, 1.0, 1.5, 2.0]
  for target_t in test_times {
    let idx = (target_t / dt).to_int()
    let y_sim = sim.outputs[idx]
    let y_analytic = @analysis.step_response_first_order(2.0, 0.5, target_t).unwrap()
    let diff = (y_sim - y_analytic).abs()
    assert_true(diff < 1.0e-6)
  }

  let max_err = error_against_first_order_analytic(2.0, 0.5, sim).unwrap()
  assert_true(max_err < 1.0e-6)
}

test "RK4 second-order step response matches analytical solution to < 1e-6" {
  // G(s) = 4 / (s^2 + 2s + 4) -> wn = 2.0, zeta = 0.5
  let tf = @tf.TransferFunction::from_coeffs([4.0], [1.0, 2.0, 4.0]).unwrap()
  let dt = 0.0005
  let sim = simulate_step(tf, 3.0, dt~).unwrap()

  let max_err2 = error_against_second_order_analytic(2.0, 0.5, sim).unwrap()
  assert_true(max_err2 < 1.0e-6)
}

test "Test 8: Closed-loop PID simulation matches second-order analytical step response to < 1e-6" {
  // Plant: G(s) = 1.0 / (s^2 + s) = 1.0 / (s(s + 1))
  // Controller: Proportional gain Kp = 4.0 (Ki = 0, Kd = 0)
  // Closed-loop transfer function T(s) = 4 / (s^2 + s + 4)
  // Standard second-order form: wn^2 = 4 -> wn = 2.0, 2*zeta*wn = 1 -> zeta = 0.25
  let plant = @tf.TransferFunction::from_coeffs([1.0], [1.0, 1.0, 0.0]).unwrap()
  let controller = @tf.TransferFunction::from_coeffs([4.0], [1.0]).unwrap()
  let closed_loop = (controller * plant).unity_feedback()

  let dt = 0.0005
  let sim = simulate_step(closed_loop, 4.0, dt~).unwrap()

  // Compare RK4 closed-loop simulation with analytic second order response
  let wn = 2.0
  let zeta = 0.25
  let check_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
  for t in check_points {
    let idx = (t / dt).to_int()
    let y_rk4 = sim.outputs[idx]
    let y_exact = @analysis.step_response_second_order(wn, zeta, t).unwrap()
    let err = (y_rk4 - y_exact).abs()
    assert_true(err < 1.0e-6)
  }
}
