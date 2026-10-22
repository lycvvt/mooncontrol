test "complex addition and subtraction" {
  let c1 = Complex::new(3.0, 4.0)
  let c2 = Complex::new(1.0, -2.0)
  let sum = c1 + c2
  assert_eq(sum.re, 4.0)
  assert_eq(sum.im, 2.0)

  let diff = c1 - c2
  assert_eq(diff.re, 2.0)
  assert_eq(diff.im, 6.0)
}

test "complex multiplication and division" {
  // (3 + 4i)(1 - 2i) = 3 - 6i + 4i - 8i^2 = 11 - 2i
  let c1 = Complex::new(3.0, 4.0)
  let c2 = Complex::new(1.0, -2.0)
  let prod = c1 * c2
  assert_eq(prod.re, 11.0)
  assert_eq(prod.im, -2.0)

  // (11 - 2i) / (1 - 2i) = (3 + 4i)
  let quot = prod / c2
  assert_true(quot.approx_eq(c1))
}

test "complex magnitude conjugate and arg" {
  let c = Complex::new(3.0, 4.0)
  assert_eq(c.abs(), 5.0)
  assert_eq(c.abs_sqr(), 25.0)

  let conj = c.conj()
  assert_eq(conj.re, 3.0)
  assert_eq(conj.im, -4.0)

  let c_real = Complex::new(2.0, 0.0)
  assert_true(c_real.is_real())
}

test "complex polar and exp" {
  let r = 2.0
  let theta = @math.PI / 3.0
  let c = Complex::from_polar(r, theta)
  assert_true((c.abs() - 2.0).abs() < 1.0e-9)
  assert_true((c.arg() - theta).abs() < 1.0e-9)

  // e^(i * pi) = -1
  let i_pi = Complex::new(0.0, @math.PI)
  let exp_res = i_pi.exp()
  assert_true(exp_res.approx_eq(Complex::new(-1.0, 0.0), epsilon=1.0e-9))
}
