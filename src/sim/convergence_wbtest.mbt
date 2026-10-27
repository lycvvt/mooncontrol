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

  // Verify that empirical convergence order is close to 4.0 (between 3.9 and 4.1)
  assert_true(study.mean_order > 3.9 && study.mean_order < 4.1)

  // Verify that a posteriori Richardson error certificate tracks actual error within 5%
  let rec = study.records[2] // dt = 0.01
  let rel_diff = (rec.richardson_error_est - rec.max_error).abs() / rec.max_error
  assert_true(rel_diff < 0.05)
}

test "Exact 2D complex-plane RK4 stability boundary and closed-loop sampled-data counterexamples" {
  // 1. On negative real axis (zeta = 1.0, pole = -1 + 0i), r_crit ≈ 2.78529356
  let real_pole = @num.Complex::new(-1.0, 0.0)
  let r_real = rk4_directional_stability_radius(real_pole)
  assert_true((r_real - 2.78529356).abs() < 1.0e-5)

  // 2. On imaginary axis (zeta = 0.0, pole = 0 + 1i), r_crit = 2*sqrt(2) ≈ 2.828427
  let imag_pole = @num.Complex::new(0.0, 1.0)
  let r_imag = rk4_directional_stability_radius(imag_pole)
  assert_true((r_imag - 2.0 * (2.0 : Double).sqrt()).abs() < 1.0e-6)

  // 3. Counterexample A (Open-loop unstable plant):
  // Plant B: Unstable Inverted Pendulum G(s) = 1 / (s^2 - 9.8)
  // Has RHP pole +3.1305 -> open-loop unstable for ALL dt > 0, so max_stable_dt == 0.0
  // (Naive 2.785 / |lambda| falsely gave 0.8896s!)
  let tf_b = @tf.TransferFunction::from_coeffs([1.0], [1.0, 0.0, -9.8]).unwrap()
  let bound_b = check_system_boundaries(tf_b)
  assert_false(bound_b.is_open_loop_stable)
  assert_true(bound_b.max_stable_dt == 0.0)

  // 4. Counterexample B (Sampled-data closed-loop PID vs open-loop plant bound):
  // Plant: G(s) = 1 / (s + 1) (open-loop pole = -1, open-loop max_stable_dt = 2.7853s)
  // Controller: Kp = 8.0, Ki = 4.0, Kd = 0.2 (continuous closed-loop is strictly Hurwitz stable!)
  // At dt = 0.5s (which is well below open-loop 2.7853s), the sampled-data ZOH + discrete PID
  // transition matrix has spectral radius rho(Phi_CL(0.5)) > 1.0 (diverges!),
  // whereas at dt = 0.01s, rho(Phi_CL(0.01)) < 1.0 (asymptotically stable!).
  // Furthermore, if unfiltered Kd = 1.5 > 1/(c*b), even dt = 0.01s has rho > 1.0 due to the ZOH
  // finite-difference derivative algebraic loop gain (c*b*Kd > 1)!
  let plant_1st = @tf.TransferFunction::from_coeffs([1.0], [1.0, 1.0]).unwrap()
  let rho_unstable_dt = closed_loop_pid_spectral_radius(plant_1st, 8.0, 4.0, 0.2, 0.5).unwrap()
  let rho_stable_dt = closed_loop_pid_spectral_radius(plant_1st, 8.0, 4.0, 0.2, 0.01).unwrap()
  let rho_unstable_kd = closed_loop_pid_spectral_radius(plant_1st, 8.0, 4.0, 1.5, 0.01).unwrap()
  assert_true(rho_unstable_dt > 1.0)
  assert_true(rho_stable_dt < 1.0)
  assert_true(rho_unstable_kd > 1.0)
}

test "Adaptive step-size RK45/Richardson solver enforces error tolerance across fast and slow plants" {
  // Fast 1st-order system with pole at -50 rad/s and gain K = 100.0:
  // G(s) = 100 / (0.02s + 1)
  // Note: A fixed step of dt = 0.01 would give |lambda|*dt = 0.5 and error ~ 0.024 >> 1e-7.
  // The adaptive solver automatically shrinks dt during the fast transient to meet tolerance!
  let fast_tf = @tf.TransferFunction::from_coeffs([100.0], [0.02, 1.0]).unwrap()
  let adap = simulate_step_adaptive(
    fast_tf,
    0.5,
    atol=1.0e-9,
    rtol=1.0e-8,
    dt_init=0.02,
  ).unwrap()

  assert_true(adap.accepted_steps > 10)
  for i = 0; i < adap.times.length(); i = i + 1 {
    let t = adap.times[i]
    let y_sim = adap.outputs[i]
    let y_exact = @analysis.step_response_first_order(100.0, 0.02, t).unwrap()
    let err = (y_sim - y_exact).abs()
    assert_true(err < 1.0e-5)
  }
}

test "Nonlinear closed-loop simulation: actuator saturation, anti-windup, transport delay, and disturbance" {
  // Plant: G(s) = 10 / (s^2 + 2s + 10)
  let plant = @tf.TransferFunction::from_coeffs([10.0], [1.0, 2.0, 10.0]).unwrap()

  // Compare saturated loop WITHOUT back-calculation anti-windup (kb = 0)
  // vs WITH back-calculation anti-windup (kb = 5.0), under actuator saturation [-1.5, 1.5],
  // pure transport delay tau_d = 0.02s, and external disturbance at t = 2.5s
  let sim_no_aw = simulate_closed_loop_nonlinear(
    plant,
    2.5,
    4.0,
    0.4,
    5.0,
    setpoint=1.0,
    dt=0.005,
    u_min=-1.5,
    u_max=1.5,
    anti_windup_kb=0.0,
    transport_delay=0.02,
    disturbance_time=2.5,
    disturbance_magnitude=-0.4,
  ).unwrap()

  let sim_with_aw = simulate_closed_loop_nonlinear(
    plant,
    2.5,
    4.0,
    0.4,
    5.0,
    setpoint=1.0,
    dt=0.005,
    u_min=-1.5,
    u_max=1.5,
    anti_windup_kb=5.0,
    transport_delay=0.02,
    disturbance_time=2.5,
    disturbance_magnitude=-0.4,
  ).unwrap()

  let metrics_no_aw = measure_trajectory_metrics(sim_no_aw, setpoint=1.0)
  let metrics_with_aw = measure_trajectory_metrics(sim_with_aw, setpoint=1.0)

  // Back-calculation anti-windup strictly reduces overshoot and IAE under saturation
  assert_true(metrics_with_aw.overshoot_percent < metrics_no_aw.overshoot_percent)
  assert_true(metrics_with_aw.iae < metrics_no_aw.iae)
  assert_true(metrics_with_aw.steady_state_error < 0.02)
}
