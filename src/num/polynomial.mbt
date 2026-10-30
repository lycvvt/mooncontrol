///|
/// Polynomial representation with coefficients in descending power order:
/// `coeffs[0]*s^n + coeffs[1]*s^(n-1) + ... + coeffs[n]`
pub(all) struct Polynomial {
  coeffs : Array[Double]
} derive(Eq)

pub extend Polynomial with Eq::{equal, not_equal}
pub extend Polynomial with Show::{output}
pub extend Polynomial with Add::{add}
pub extend Polynomial with Sub::{sub}
pub extend Polynomial with Mul::{mul}

///|
/// Create a new Polynomial, trimming any superfluous leading zeros.
pub fn Polynomial::new(coeffs : Array[Double]) -> Polynomial {
  if coeffs.length() == 0 {
    return { coeffs: [0.0] }
  }
  let mut start = 0
  while start < coeffs.length() - 1 && coeffs[start].abs() < 1.0e-12 {
    start = start + 1
  }
  let trimmed = Array(capacity=coeffs.length() - start)
  for i = start; i < coeffs.length(); i = i + 1 {
    trimmed.push(coeffs[i])
  }
  { coeffs: trimmed }
}

///|
/// Polynomial degree (e.g. [3.0, 2.0] is degree 1).
pub fn Polynomial::degree(self : Polynomial) -> Int {
  self.coeffs.length() - 1
}

///|
/// Check if polynomial is approximately zero.
pub fn Polynomial::is_zero(self : Polynomial, epsilon~ : Double = 1.0e-12) -> Bool {
  for c in self.coeffs {
    if c.abs() > epsilon {
      return false
    }
  }
  true
}

///|
/// Evaluate polynomial at real value x using Horner's method.
pub fn Polynomial::eval(self : Polynomial, x : Double) -> Double {
  let mut acc = self.coeffs[0]
  for i = 1; i < self.coeffs.length(); i = i + 1 {
    acc = acc * x + self.coeffs[i]
  }
  acc
}

///|
/// Evaluate polynomial at complex value z using Horner's method.
pub fn Polynomial::eval_complex(self : Polynomial, z : Complex) -> Complex {
  let mut acc = Complex::from_real(self.coeffs[0])
  for i = 1; i < self.coeffs.length(); i = i + 1 {
    acc = acc * z + Complex::from_real(self.coeffs[i])
  }
  acc
}

///|
/// Multiply polynomial by scalar.
pub fn Polynomial::scale(self : Polynomial, factor : Double) -> Polynomial {
  let res = Array(capacity=self.coeffs.length())
  for c in self.coeffs {
    res.push(c * factor)
  }
  Polynomial::new(res)
}

///|
/// Addition of two polynomials.
pub impl Add for Polynomial with add(self, other) {
  let len1 = self.coeffs.length()
  let len2 = other.coeffs.length()
  let max_len = if len1 > len2 { len1 } else { len2 }
  let res = Array::make(max_len, 0.0)
  for i = 0; i < max_len; i = i + 1 {
    let idx1 = len1 - 1 - i
    let idx2 = len2 - 1 - i
    let c1 = if idx1 >= 0 { self.coeffs[idx1] } else { 0.0 }
    let c2 = if idx2 >= 0 { other.coeffs[idx2] } else { 0.0 }
    res[max_len - 1 - i] = c1 + c2
  }
  Polynomial::new(res)
}

///|
/// Subtraction of two polynomials.
pub impl Sub for Polynomial with sub(self, other) {
  let len1 = self.coeffs.length()
  let len2 = other.coeffs.length()
  let max_len = if len1 > len2 { len1 } else { len2 }
  let res = Array::make(max_len, 0.0)
  for i = 0; i < max_len; i = i + 1 {
    let idx1 = len1 - 1 - i
    let idx2 = len2 - 1 - i
    let c1 = if idx1 >= 0 { self.coeffs[idx1] } else { 0.0 }
    let c2 = if idx2 >= 0 { other.coeffs[idx2] } else { 0.0 }
    res[max_len - 1 - i] = c1 - c2
  }
  Polynomial::new(res)
}

///|
/// Multiplication of two polynomials (convolution).
pub impl Mul for Polynomial with mul(self, other) {
  let len1 = self.coeffs.length()
  let len2 = other.coeffs.length()
  let out_len = len1 + len2 - 1
  let res = Array::make(out_len, 0.0)
  for i = 0; i < len1; i = i + 1 {
    for j = 0; j < len2; j = j + 1 {
      res[i + j] = res[i + j] + self.coeffs[i] * other.coeffs[j]
    }
  }
  Polynomial::new(res)
}

///|
/// Formal derivative of polynomial.
pub fn Polynomial::derivative(self : Polynomial) -> Polynomial {
  let deg = self.degree()
  if deg <= 0 {
    return Polynomial::new([0.0])
  }
  let res = Array(capacity=deg)
  for i = 0; i < deg; i = i + 1 {
    let power = (deg - i).to_double()
    res.push(self.coeffs[i] * power)
  }
  Polynomial::new(res)
}

///|
/// Find all roots of the polynomial.
/// Degree 1: analytic
/// Degree 2: quadratic formula
/// Degree >= 3: Durand-Kerner simultaneous iteration
pub fn Polynomial::roots(self : Polynomial) -> Array[Complex] {
  let deg = self.degree()
  if deg <= 0 {
    return []
  }
  if deg == 1 {
    let a1 = self.coeffs[0]
    let a0 = self.coeffs[1]
    return [Complex::from_real(-a0 / a1)]
  }
  if deg == 2 {
    let a = self.coeffs[0]
    let b = self.coeffs[1]
    let c = self.coeffs[2]
    let disc = b * b - 4.0 * a * c
    if disc >= 0.0 {
      let sqrt_d = disc.sqrt()
      let r1 = (-b + sqrt_d) / (2.0 * a)
      let r2 = (-b - sqrt_d) / (2.0 * a)
      return [Complex::from_real(r1), Complex::from_real(r2)]
    } else {
      let sqrt_neg_d = (-disc).sqrt()
      let re = -b / (2.0 * a)
      let im = sqrt_neg_d / (2.0 * a)
      return [Complex::new(re, im), Complex::new(re, -im)]
    }
  }

  // Deflate exact roots at s = 0 (trailing zero coefficients) first
  let mut zero_roots = 0
  let mut last_idx = self.coeffs.length() - 1
  while last_idx > 0 && self.coeffs[last_idx].abs() < 1.0e-14 {
    zero_roots = zero_roots + 1
    last_idx = last_idx - 1
  }
  if zero_roots > 0 {
    let deflated_coeffs = Array(capacity=last_idx + 1)
    for i = 0; i <= last_idx; i = i + 1 {
      deflated_coeffs.push(self.coeffs[i])
    }
    let deflated_poly = Polynomial::new(deflated_coeffs)
    let all_roots = deflated_poly.roots()
    for z = 0; z < zero_roots; z = z + 1 {
      all_roots.push(Complex::zero())
    }
    return all_roots
  }

  // Higher degree: Durand-Kerner method with geometric root-radius scaling
  // Normalize polynomial so leading coefficient is 1.0
  let lead = self.coeffs[0]
  let monic_coeffs = Array(capacity=self.coeffs.length())
  for c in self.coeffs {
    monic_coeffs.push(c / lead)
  }
  let monic = Polynomial::new(monic_coeffs)

  let n = deg
  let const_mag = monic.coeffs[n].abs()
  let radius = if const_mag > 1.0e-12 {
    let r_geom = @math.pow(const_mag, 1.0 / n.to_double())
    if r_geom > 0.5 { r_geom } else { 0.5 }
  } else {
    1.0
  }

  let current_roots = Array::make(n, Complex::zero())
  // Initial points distributed around circle of radius `radius` with asymmetric phase offset
  let two_pi = 8.0 * @math.atan(1.0)
  for k = 0; k < n; k = k + 1 {
    let angle = two_pi * k.to_double() / n.to_double() + 0.37
    current_roots[k] = Complex::from_polar(radius, angle)
  }

  let max_iter = 400
  for _iter = 0; _iter < max_iter; _iter = _iter + 1 {
    let mut max_diff = 0.0
    for k = 0; k < n; k = k + 1 {
      let zk = current_roots[k]
      let p_val = monic.eval_complex(zk)
      let mut denom = Complex::one()
      for j = 0; j < n; j = j + 1 {
        if j != k {
          denom = denom * (zk - current_roots[j])
        }
      }
      let delta = if denom.abs_sqr() > 1.0e-24 {
        p_val / denom
      } else {
        Complex::zero()
      }
      let new_zk = zk - delta
      let diff = delta.abs()
      if diff > max_diff {
        max_diff = diff
      }
      current_roots[k] = new_zk
    }
    if max_diff < 1.0e-13 {
      break
    }
  }

  // Clean up tiny imaginary parts
  for k = 0; k < n; k = k + 1 {
    let root = current_roots[k]
    if root.im.abs() < 1.0e-10 {
      current_roots[k] = Complex::from_real(root.re)
    }
  }

  current_roots
}

///|
pub fn Polynomial::to_string(self : Polynomial) -> String {
  let deg = self.degree()
  if deg == 0 {
    return "\{self.coeffs[0]}"
  }
  let buf = StringBuilder()
  for i = 0; i <= deg; i = i + 1 {
    let c = self.coeffs[i]
    if c == 0.0 && deg > 0 {
      continue
    }
    let power = deg - i
    if buf.is_empty() {
      if power == 0 {
        buf.write_string("\{c}")
      } else if power == 1 {
        if c == 1.0 {
          buf.write_string("s")
        } else if c == -1.0 {
          buf.write_string("-s")
        } else {
          buf.write_string("\{c}s")
        }
      } else if c == 1.0 {
        buf.write_string("s^\{power}")
      } else if c == -1.0 {
        buf.write_string("-s^\{power}")
      } else {
        buf.write_string("\{c}s^\{power}")
      }
    } else {
      let sign_str = if c > 0.0 { " + " } else { " - " }
      buf.write_string(sign_str)
      let abs_c = c.abs()
      if power == 0 {
        buf.write_string("\{abs_c}")
      } else if power == 1 {
        if abs_c == 1.0 {
          buf.write_string("s")
        } else {
          buf.write_string("\{abs_c}s")
        }
      } else if abs_c == 1.0 {
        buf.write_string("s^\{power}")
      } else {
        buf.write_string("\{abs_c}s^\{power}")
      }
    }
  }
  if buf.is_empty() {
    "0"
  } else {
    buf.to_string()
  }
}

///|
pub impl Show for Polynomial with output(self, logger) {
  logger.write_string(self.to_string())
}
