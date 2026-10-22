///|
/// Complex number representation with double-precision real and imaginary parts.
pub(all) struct Complex {
  re : Double
  im : Double
} derive(Eq)

pub extend Complex with Eq::{equal, not_equal}
pub extend Complex with Show::{output}
pub extend Complex with Add::{add}
pub extend Complex with Sub::{sub}
pub extend Complex with Neg::{neg}
pub extend Complex with Mul::{mul}
pub extend Complex with Div::{div}

///|
pub fn Complex::new(re : Double, im : Double) -> Complex {
  { re, im }
}

///|
pub fn Complex::from_real(re : Double) -> Complex {
  { re, im: 0.0 }
}

///|
pub fn Complex::zero() -> Complex {
  { re: 0.0, im: 0.0 }
}

///|
pub fn Complex::one() -> Complex {
  { re: 1.0, im: 0.0 }
}

///|
pub fn Complex::i() -> Complex {
  { re: 0.0, im: 1.0 }
}

///|
pub fn Complex::from_polar(r : Double, theta : Double) -> Complex {
  { re: r * @math.cos(theta), im: r * @math.sin(theta) }
}

///|
pub impl Add for Complex with add(self, other) {
  { re: self.re + other.re, im: self.im + other.im }
}

///|
pub impl Sub for Complex with sub(self, other) {
  { re: self.re - other.re, im: self.im - other.im }
}

///|
pub impl Neg for Complex with neg(self) {
  { re: -self.re, im: -self.im }
}

///|
pub impl Mul for Complex with mul(self, other) {
  {
    re: self.re * other.re - self.im * other.im,
    im: self.re * other.im + self.im * other.re,
  }
}

///|
pub impl Div for Complex with div(self, other) {
  let denom = other.re * other.re + other.im * other.im
  {
    re: (self.re * other.re + self.im * other.im) / denom,
    im: (self.im * other.re - self.re * other.im) / denom,
  }
}

///|
pub fn Complex::add_real(self : Complex, r : Double) -> Complex {
  { re: self.re + r, im: self.im }
}

///|
pub fn Complex::sub_real(self : Complex, r : Double) -> Complex {
  { re: self.re - r, im: self.im }
}

///|
pub fn Complex::mul_real(self : Complex, r : Double) -> Complex {
  { re: self.re * r, im: self.im * r }
}

///|
pub fn Complex::div_real(self : Complex, r : Double) -> Complex {
  { re: self.re / r, im: self.im / r }
}

///|
pub fn Complex::abs(self : Complex) -> Double {
  @math.hypot(self.re, self.im)
}

///|
pub fn Complex::abs_sqr(self : Complex) -> Double {
  self.re * self.re + self.im * self.im
}

///|
pub fn Complex::arg(self : Complex) -> Double {
  @math.atan2(self.im, self.re)
}

///|
pub fn Complex::conj(self : Complex) -> Complex {
  { re: self.re, im: -self.im }
}

///|
pub fn Complex::exp(self : Complex) -> Complex {
  let r = @math.exp(self.re)
  { re: r * @math.cos(self.im), im: r * @math.sin(self.im) }
}

///|
pub fn Complex::sqrt(self : Complex) -> Complex {
  let r = self.abs()
  let gamma = ((r + self.re) / 2.0).sqrt()
  let delta = if self.im < 0.0 {
    -(((r - self.re) / 2.0).sqrt())
  } else {
    ((r - self.re) / 2.0).sqrt()
  }
  { re: gamma, im: delta }
}

///|
pub fn Complex::approx_eq(
  self : Complex,
  other : Complex,
  epsilon~ : Double = 0.000000001
) -> Bool {
  (self.re - other.re).abs() <= epsilon && (self.im - other.im).abs() <= epsilon
}

///|
pub fn Complex::is_real(self : Complex, epsilon~ : Double = 0.000000001) -> Bool {
  self.im.abs() <= epsilon
}

///|
pub fn Complex::to_string(self : Complex) -> String {
  if self.im == 0.0 {
    "\{self.re}"
  } else if self.re == 0.0 {
    "\{self.im}i"
  } else if self.im > 0.0 {
    "\{self.re}+\{self.im}i"
  } else {
    "\{self.re}\{self.im}i"
  }
}

///|
pub impl Show for Complex with output(self, logger) {
  logger.write_string(self.to_string())
}
