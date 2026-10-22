test "pid proportional and integral" {
  let pid = PIDController::new(2.0, 1.0, 0.0)
  // Step error of 1.0 at dt = 0.1
  // Step 1: P = 2.0 * 1.0 = 2.0; I = 1.0 * (1.0 * 0.1) = 0.1 -> u = 2.1
  let u1 = pid.update(1.0, 0.0, 0.1)
  assert_true((u1 - 2.1).abs() < 1.0e-9)

  // Step 2: P = 2.0; I = 1.0 * 0.2 = 0.2 -> u = 2.2
  let u2 = pid.update(1.0, 0.0, 0.1)
  assert_true((u2 - 2.2).abs() < 1.0e-9)
}

test "pid anti-windup clamping" {
  // Integral limit = 0.5
  let pid = PIDController::new(0.0, 1.0, 0.0, integral_limit=0.5)
  // 10 steps of error = 1.0, dt = 0.1 -> integral without limit would be 1.0, but clamped to 0.5
  for i = 0; i < 10; i = i + 1 {
    let _ = pid.update(1.0, 0.0, 0.1)
  }
  assert_true((pid.integral - 0.5).abs() < 1.0e-9)
}

test "pid derivative action" {
  let pid = PIDController::new(0.0, 0.0, 1.0)
  // First update initializes prev_error, derivative is 0
  let u1 = pid.update(1.0, 0.0, 0.1)
  assert_eq(u1, 0.0)

  // Second update: error changes from 1.0 to 1.5 -> de/dt = 0.5 / 0.1 = 5.0 -> D = 5.0
  let u2 = pid.update(1.5, 0.0, 0.1)
  assert_true((u2 - 5.0).abs() < 1.0e-9)
}

test "pid to transfer function" {
  let pid = PIDController::new(3.0, 2.0, 1.0)
  let tf = pid.to_transfer_function().unwrap()
  // C(s) = (s^2 + 3s + 2) / s
  assert_eq(tf.num.coeffs, [1.0, 3.0, 2.0])
  assert_eq(tf.den.coeffs, [1.0, 0.0])
}
