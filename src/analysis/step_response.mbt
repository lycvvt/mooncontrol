///|
/// Analytic unit step response of standard first-order system:
/// G(s) = K / (T*s + 1)
/// y(t) = K * (1 - e^(-t / T)) for t >= 0
pub fn step_response_first_order(
  k : Double,
  t_const : Double,
  t : Double
) -> Result[Double, String] {
  if t_const <= 0.0 {
    return Err("Time constant T must be strictly positive")
  }
  if t < 0.0 {
    return Ok(0.0)
  }
  let exponent = -t / t_const
  let y = k * (1.0 - @math.exp(exponent))
  Ok(y)
}

///|
/// Analytic unit step response of standard second-order system:
/// G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
/// Covers underdamped (0 < zeta < 1), critically damped (zeta = 1), and overdamped (zeta > 1).
pub fn step_response_second_order(
  omega_n : Double,
  zeta : Double,
  t : Double
) -> Result[Double, String] {
  if omega_n <= 0.0 {
    return Err("Natural frequency omega_n must be strictly positive")
  }
  if zeta < 0.0 {
    return Err("Damping ratio zeta must be non-negative for stable system")
  }
  if t < 0.0 {
    return Ok(0.0)
  }
  if t == 0.0 {
    return Ok(0.0)
  }

  // Underdamped case: 0 <= zeta < 1
  if zeta < 0.999999999 {
    let sqrt_1_minus_zeta2 = (1.0 - zeta * zeta).sqrt()
    let omega_d = omega_n * sqrt_1_minus_zeta2
    let phi = @math.acos(zeta)
    let decay = @math.exp(-zeta * omega_n * t)
    let oscillation = @math.sin(omega_d * t + phi)
    let y = 1.0 - decay / sqrt_1_minus_zeta2 * oscillation
    return Ok(y)
  }

  // Critically damped case: zeta ≈ 1
  if zeta <= 1.000000001 {
    let decay = @math.exp(-omega_n * t)
    let y = 1.0 - decay * (1.0 + omega_n * t)
    return Ok(y)
  }

  // Overdamped case: zeta > 1
  let sqrt_zeta2_minus_1 = (zeta * zeta - 1.0).sqrt()
  let s1 = -omega_n * (zeta - sqrt_zeta2_minus_1)
  let s2 = -omega_n * (zeta + sqrt_zeta2_minus_1)
  let denom = s1 - s2
  let term1 = s2 / denom * @math.exp(s1 * t)
  let term2 = s1 / denom * @math.exp(s2 * t)
  let y = 1.0 + term1 - term2
  Ok(y)
}
