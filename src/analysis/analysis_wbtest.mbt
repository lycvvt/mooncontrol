test "Test 3: First-order analytic step response matches K*(1-e^(-t/T)) to 1e-9" {
  // Plant: G(s) = 2.0 / (0.5s + 1.0) -> K = 2.0, T = 0.5
  let k = 2.0
  let t_const = 0.5

  // At t = 0: y(0) = 0
  let y0 = step_response_first_order(k, t_const, 0.0).unwrap()
  assert_true(y0.abs() < 1.0e-9)

  // At t = T: 1 - e^(-1) ≈ 0.6321205588285577 -> K*(1 - e^(-1))
  let y_t = step_response_first_order(k, t_const, t_const).unwrap()
  let expected_t = k * (1.0 - @math.exp(-1.0))
  assert_true((y_t - expected_t).abs() < 1.0e-9)

  // At t = 2T: K*(1 - e^(-2))
  let y_2t = step_response_first_order(k, t_const, 2.0 * t_const).unwrap()
  let expected_2t = k * (1.0 - @math.exp(-2.0))
  assert_true((y_2t - expected_2t).abs() < 1.0e-9)

  // At t = 3T: K*(1 - e^(-3))
  let y_3t = step_response_first_order(k, t_const, 3.0 * t_const).unwrap()
  let expected_3t = k * (1.0 - @math.exp(-3.0))
  assert_true((y_3t - expected_3t).abs() < 1.0e-9)
}

test "Test 4: Second-order analytic step response for underdamped, critical, and overdamped" {
  let omega_n = 2.0

  // 1. Underdamped (zeta = 0.5)
  // G(s) = 4 / (s^2 + 2s + 4)
  let zeta_under = 0.5
  let t_under = 1.2
  let y_under = step_response_second_order(omega_n, zeta_under, t_under).unwrap()
  let wd = omega_n * (1.0 - zeta_under * zeta_under).sqrt()
  let phi = @math.acos(zeta_under)
  let expected_under = 1.0 - @math.exp(-zeta_under * omega_n * t_under) / (1.0 - zeta_under * zeta_under).sqrt() * @math.sin(wd * t_under + phi)
  assert_true((y_under - expected_under).abs() < 1.0e-9)

  // 2. Critically damped (zeta = 1.0)
  // G(s) = 4 / (s^2 + 4s + 4)
  let zeta_crit = 1.0
  let t_crit = 0.8
  let y_crit = step_response_second_order(omega_n, zeta_crit, t_crit).unwrap()
  let expected_crit = 1.0 - @math.exp(-omega_n * t_crit) * (1.0 + omega_n * t_crit)
  assert_true((y_crit - expected_crit).abs() < 1.0e-9)

  // 3. Overdamped (zeta = 2.0)
  // G(s) = 4 / (s^2 + 8s + 4)
  let zeta_over = 2.0
  let t_over = 1.5
  let y_over = step_response_second_order(omega_n, zeta_over, t_over).unwrap()
  let sq_zeta = (zeta_over * zeta_over - 1.0).sqrt()
  let p1 = -omega_n * (zeta_over - sq_zeta)
  let p2 = -omega_n * (zeta_over + sq_zeta)
  let expected_over = 1.0 + p2 / (p1 - p2) * @math.exp(p1 * t_over) - p1 / (p1 - p2) * @math.exp(p2 * t_over)
  assert_true((y_over - expected_over).abs() < 1.0e-9)

  // At t = 0, all responses must be exactly 0
  assert_true(step_response_second_order(omega_n, zeta_under, 0.0).unwrap().abs() < 1.0e-9)
  assert_true(step_response_second_order(omega_n, zeta_crit, 0.0).unwrap().abs() < 1.0e-9)
  assert_true(step_response_second_order(omega_n, zeta_over, 0.0).unwrap().abs() < 1.0e-9)
}

test "Test 5: Performance metrics Mp, tp, ts vs theoretical formulas < 1e-6" {
  // Standard test case: omega_n = 5.0, zeta = 0.6
  let wn = 5.0
  let zeta = 0.6
  let metrics = calculate_second_order_metrics(wn, zeta).unwrap()

  // Mp = exp(-pi * 0.6 / sqrt(1 - 0.36)) = exp(-0.6*pi / 0.8) = exp(-0.75*pi) ≈ 0.09478
  let expected_mp = @math.exp(-0.75 * @math.PI)
  assert_true((metrics.overshoot_fraction - expected_mp).abs() < 1.0e-6)
  assert_true((metrics.overshoot_percent - expected_mp * 100.0).abs() < 1.0e-6)

  // tp = pi / (wn * sqrt(1 - 0.36)) = pi / (5 * 0.8) = pi / 4.0 ≈ 0.785398
  let expected_tp = @math.PI / 4.0
  assert_true((metrics.peak_time - expected_tp).abs() < 1.0e-6)

  // ts_2pct = 4 / (0.6 * 5) = 4 / 3 ≈ 1.333333
  let expected_ts = 4.0 / 3.0
  assert_true((metrics.settling_time_2pct - expected_ts).abs() < 1.0e-6)
}

test "Test 6: Steady-state error for Type 0 and Type 1 systems vs final value theorem" {
  // Type 0 system: G(s) = 5 / (s + 2) -> Kp = 5 / 2 = 2.5
  let g_type0 = @tf.TransferFunction::from_coeffs([5.0], [1.0, 2.0]).unwrap()
  // Step input: e_ss = 1 / (1 + 2.5) = 1 / 3.5 ≈ 0.2857142857
  let ess_step0 = steady_state_error(g_type0, Step).unwrap()
  assert_true((ess_step0 - 1.0 / 3.5).abs() < 1.0e-9)

  // Ramp input: e_ss = infinity for Type 0
  let ess_ramp0 = steady_state_error(g_type0, Ramp).unwrap()
  assert_eq(ess_ramp0, @double.infinity)

  // Type 1 system: G(s) = 10 / (s(s + 2)) = 10 / (s^2 + 2s) -> Kv = 10 / 2 = 5.0
  let g_type1 = @tf.TransferFunction::from_coeffs([10.0], [1.0, 2.0, 0.0]).unwrap()
  // Step input: e_ss = 0 for Type 1
  let ess_step1 = steady_state_error(g_type1, Step).unwrap()
  assert_true(ess_step1.abs() < 1.0e-9)

  // Ramp input: e_ss = 1 / Kv = 1 / 5 = 0.2
  let ess_ramp1 = steady_state_error(g_type1, Ramp).unwrap()
  assert_true((ess_ramp1 - 0.2).abs() < 1.0e-9)
}

test "Test 7: Routh-Hurwitz stability criterion agrees with pole real-part criterion" {
  // Case 1: Stable system s^3 + 6s^2 + 11s + 6 = (s+1)(s+2)(s+3)
  let p_stable = @num.Polynomial::new([1.0, 6.0, 11.0, 6.0])
  let tf_stable = @tf.TransferFunction::from_coeffs([1.0], [1.0, 6.0, 11.0, 6.0]).unwrap()
  let rh_stable = routh_hurwitz(p_stable).unwrap()
  assert_true(rh_stable.is_stable)
  assert_true(is_stable_by_poles(tf_stable))
  assert_eq(rh_stable.is_stable, is_stable_by_poles(tf_stable))

  // Case 2: Unstable system s^3 + s^2 + 2s + 24 (has right half plane roots)
  let p_unstable = @num.Polynomial::new([1.0, 1.0, 2.0, 24.0])
  let tf_unstable = @tf.TransferFunction::from_coeffs([1.0], [1.0, 1.0, 2.0, 24.0]).unwrap()
  let rh_unstable = routh_hurwitz(p_unstable).unwrap()
  assert_true(!rh_unstable.is_stable)
  assert_true(!is_stable_by_poles(tf_unstable))
  assert_eq(rh_unstable.is_stable, is_stable_by_poles(tf_unstable))
  assert_true(rh_unstable.sign_changes > 0)

  // Case 3: Stable 4th order system s^4 + 3s^3 + 3s^2 + 2s + 1
  // Routh array:
  // s^4: 1   3   1
  // s^3: 3   2   0
  // s^2: (9-2)/3 = 7/3, 1
  // s^1: (14/3 - 3)/(7/3) = 5/7
  // s^0: 1
  // All positive -> stable!
  let p4 = @num.Polynomial::new([1.0, 3.0, 3.0, 2.0, 1.0])
  let tf4 = @tf.TransferFunction::from_coeffs([1.0], [1.0, 3.0, 3.0, 2.0, 1.0]).unwrap()
  let rh4 = routh_hurwitz(p4).unwrap()
  assert_true(rh4.is_stable)
  assert_true(is_stable_by_poles(tf4))
  assert_eq(rh4.is_stable, is_stable_by_poles(tf4))
}
