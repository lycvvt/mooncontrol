test "RK4 4th-order convergence on 2nd-order underdamped DC motor plant" {
  // Plant A: G(s) = 10 / (s^2 + 2s + 10)
  // wn = sqrt(10) ≈ 3.16227766, zeta = 2 / (2 * wn) = 1 / sqrt(10) ≈ 0.316227766
  let wn = (10.0 : Double).sqrt()
  let zeta = 1.0 / wn
  let tf = @tf.TransferFunction::from_coeffs([10.0], [1.0, 2.0, 10.0]).unwrap()

  // Sweep step size dt by halving: 0.04, 0.02, 0.01, 0.005
  let dts = [0.04, 0.02, 0.01, 0.005]
  let duration = 3.0
  let study = evaluate_rk4_convergence(tf, duration, dts, fn(t) {
    @analysis.step_response_second_order(wn, zeta, t).unwrap()
  }).unwrap()

  // Assert that errors strictly decrease as dt decreases
  assert_true(study.records[0].max_error > study.records[1].max_error)
  assert_true(study.records[1].max_error > study.records[2].max_error)
  assert_true(study.records[2].max_error > study.records[3].max_error)

  // Verify that empirical convergence order is close to 4.0 (between 3.8 and 4.2)
  assert_true(study.mean_order > 3.8 && study.mean_order < 4.2)

  // At dt = 0.005, error is already down into 10^-8 range
  assert_true(study.records[3].max_error < 1.0e-7)
}

test "RK4 4th-order convergence on 1st-order system" {
  // G(s) = 2.0 / (0.5s + 1.0) -> K = 2.0, T = 0.5
  let tf = @tf.TransferFunction::from_coeffs([2.0], [0.5, 1.0]).unwrap()
  let dts = [0.08, 0.04, 0.02, 0.01]
  let duration = 2.5
  let study = evaluate_rk4_convergence(tf, duration, dts, fn(t) {
    @analysis.step_response_first_order(2.0, 0.5, t).unwrap()
  }).unwrap()

  // Mean order must be close to 4.0
  assert_true(study.mean_order > 3.8 && study.mean_order < 4.2)
  assert_true(study.records[3].max_error < 1.0e-8)
}

test "System stability boundary and RK4 step size limit analysis" {
  // Plant A: G(s) = 10 / (s^2 + 2s + 10)
  let tf_a = @tf.TransferFunction::from_coeffs([10.0], [1.0, 2.0, 10.0]).unwrap()
  let bound_a = check_system_boundaries(tf_a)

  assert_true(bound_a.is_open_loop_stable)
  // Pole magnitude is sqrt(10) ≈ 3.162
  assert_true((bound_a.max_pole_magnitude - (10.0 : Double).sqrt()).abs() < 1.0e-4)
  // max_stable_dt ≈ 2.785 / 3.162 ≈ 0.88s
  assert_true(bound_a.max_stable_dt > 0.8 && bound_a.max_stable_dt < 1.0)

  // RK4 amplification factor at safe step dt = 0.01 must be strictly <= 1.0
  let pole = tf_a.poles()[0]
  let amp_safe = rk4_amplification_factor(pole, 0.01)
  assert_true(amp_safe <= 1.0)

  // Plant B: Unstable Inverted Pendulum G(s) = 1 / (s^2 - 9.8)
  let tf_b = @tf.TransferFunction::from_coeffs([1.0], [1.0, 0.0, -9.8]).unwrap()
  let bound_b = check_system_boundaries(tf_b)
  assert_false(bound_b.is_open_loop_stable)
}
