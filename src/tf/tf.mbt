///|
/// Continuous-time Linear Time-Invariant (LTI) Transfer Function G(s) = num(s) / den(s)
pub(all) struct TransferFunction {
  num : @num.Polynomial
  den : @num.Polynomial
} derive(Eq)

pub extend TransferFunction with Eq::{equal, not_equal}
pub extend TransferFunction with Show::{output}
pub extend TransferFunction with Add::{add}
pub extend TransferFunction with Sub::{sub}
pub extend TransferFunction with Mul::{mul}
pub extend TransferFunction with Div::{div}

///|
/// Create a new TransferFunction with validation and monic normalization.
pub fn TransferFunction::new(
  num : @num.Polynomial,
  den : @num.Polynomial
) -> Result[TransferFunction, String] {
  if den.is_zero() {
    return Err("Denominator of transfer function cannot be zero polynomial")
  }
  let lead = den.coeffs[0]
  if lead.abs() < 1.0e-12 {
    return Err("Denominator leading coefficient is zero")
  }
  // Monic normalization: divide both num and den by den's leading coefficient
  let norm_num = num.scale(1.0 / lead)
  let norm_den = den.scale(1.0 / lead)
  Ok({ num: norm_num, den: norm_den })
}

///|
/// Create a TransferFunction from raw coefficient arrays.
pub fn TransferFunction::from_coeffs(
  num : Array[Double],
  den : Array[Double]
) -> Result[TransferFunction, String] {
  let p_num = @num.Polynomial::new(num)
  let p_den = @num.Polynomial::new(den)
  TransferFunction::new(p_num, p_den)
}

///|
/// Standard first-order system: G(s) = K / (T*s + 1)
pub fn TransferFunction::first_order(
  k : Double,
  t : Double
) -> Result[TransferFunction, String] {
  if t.abs() < 1.0e-12 {
    return Err("Time constant T cannot be zero")
  }
  TransferFunction::from_coeffs([k], [t, 1.0])
}

///|
/// Standard second-order system: G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
pub fn TransferFunction::second_order(
  omega_n : Double,
  zeta : Double
) -> Result[TransferFunction, String] {
  if omega_n <= 0.0 {
    return Err("Natural frequency omega_n must be positive")
  }
  let wn2 = omega_n * omega_n
  let two_zeta_wn = 2.0 * zeta * omega_n
  TransferFunction::from_coeffs([wn2], [1.0, two_zeta_wn, wn2])
}

///|
/// System order (degree of denominator polynomial).
pub fn TransferFunction::order(self : TransferFunction) -> Int {
  self.den.degree()
}

///|
/// Series connection: G1 * G2 = (num1 * num2) / (den1 * den2)
pub impl Mul for TransferFunction with mul(self, other) {
  let new_num = self.num * other.num
  let new_den = self.den * other.den
  match TransferFunction::new(new_num, new_den) {
    Ok(tf) => tf
    Err(_) => { num: new_num, den: new_den }
  }
}

///|
/// Parallel connection: G1 + G2 = (num1 * den2 + num2 * den1) / (den1 * den2)
pub impl Add for TransferFunction with add(self, other) {
  let new_num = self.num * other.den + other.num * self.den
  let new_den = self.den * other.den
  match TransferFunction::new(new_num, new_den) {
    Ok(tf) => tf
    Err(_) => { num: new_num, den: new_den }
  }
}

///|
/// Parallel subtraction: G1 - G2
pub impl Sub for TransferFunction with sub(self, other) {
  let new_num = self.num * other.den - other.num * self.den
  let new_den = self.den * other.den
  match TransferFunction::new(new_num, new_den) {
    Ok(tf) => tf
    Err(_) => { num: new_num, den: new_den }
  }
}

///|
/// Division of transfer functions: G1 / G2 = (num1 * den2) / (den1 * num2)
pub impl Div for TransferFunction with div(self, other) {
  let new_num = self.num * other.den
  let new_den = self.den * other.num
  match TransferFunction::new(new_num, new_den) {
    Ok(tf) => tf
    Err(_) => { num: new_num, den: new_den }
  }
}

///|
/// Scalar scaling: k * G(s)
pub fn TransferFunction::scale(self : TransferFunction, factor : Double) -> TransferFunction {
  let new_num = self.num.scale(factor)
  { num: new_num, den: self.den }
}

///|
/// Closed-loop feedback interconnection:
/// G_cl(s) = G / (1 - sign * G * H) = (num_G * den_H) / (den_G * den_H - sign * num_G * num_H)
/// Default sign is -1.0 for negative feedback: G / (1 + G*H)
pub fn TransferFunction::feedback(
  self : TransferFunction,
  h : TransferFunction,
  sign~ : Double = -1.0
) -> TransferFunction {
  let forward_num = self.num * h.den
  let loop_term = (self.num * h.num).scale(-sign)
  let closed_den = self.den * h.den + loop_term
  match TransferFunction::new(forward_num, closed_den) {
    Ok(tf) => tf
    Err(_) => { num: forward_num, den: closed_den }
  }
}

///|
/// Unity negative feedback closed loop: G / (1 + G)
pub fn TransferFunction::unity_feedback(self : TransferFunction) -> TransferFunction {
  let unity = TransferFunction::from_coeffs([1.0], [1.0]).unwrap()
  self.feedback(unity, sign=-1.0)
}

///|
/// System zeros (roots of numerator polynomial).
pub fn TransferFunction::zeros(self : TransferFunction) -> Array[@num.Complex] {
  self.num.roots()
}

///|
/// System poles (roots of denominator polynomial).
pub fn TransferFunction::poles(self : TransferFunction) -> Array[@num.Complex] {
  self.den.roots()
}

///|
/// Evaluate frequency response / complex transfer function at s.
pub fn TransferFunction::eval_complex(
  self : TransferFunction,
  s : @num.Complex
) -> @num.Complex {
  let num_val = self.num.eval_complex(s)
  let den_val = self.den.eval_complex(s)
  num_val / den_val
}

///|
/// DC gain (s -> 0): num(0) / den(0)
pub fn TransferFunction::dc_gain(self : TransferFunction) -> Double? {
  let den0 = self.den.eval(0.0)
  if den0.abs() < 1.0e-12 {
    None
  } else {
    Some(self.num.eval(0.0) / den0)
  }
}

///|
pub fn TransferFunction::to_string(self : TransferFunction) -> String {
  let num_str = self.num.to_string()
  let den_str = self.den.to_string()
  let max_len = if num_str.length() > den_str.length() {
    num_str.length()
  } else {
    den_str.length()
  }
  let bar = String::make(max_len + 4, '-')
  "\n  \{num_str}\n  \{bar}\n  \{den_str}\n"
}

///|
pub impl Show for TransferFunction with output(self, logger) {
  logger.write_string(self.to_string())
}
