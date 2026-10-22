test "polynomial evaluation real and complex" {
  // P(s) = 2s^2 - 3s + 1
  let p = Polynomial::new([2.0, -3.0, 1.0])
  assert_eq(p.eval(2.0), 3.0)
  assert_eq(p.eval(0.0), 1.0)
  assert_eq(p.eval(1.0), 0.0)

  // P(i) = 2(-1) - 3i + 1 = -1 - 3i
  let res_i = p.eval_complex(Complex::i())
  assert_true(res_i.approx_eq(Complex::new(-1.0, -3.0)))
}

test "polynomial addition subtraction and multiplication" {
  let p1 = Polynomial::new([1.0, 2.0, 3.0]) // s^2 + 2s + 3
  let p2 = Polynomial::new([3.0, 4.0])      // 3s + 4

  let sum = p1 + p2
  assert_eq(sum.coeffs, [1.0, 5.0, 7.0])

  let diff = p1 - p2
  assert_eq(diff.coeffs, [1.0, -1.0, -1.0])

  // (s + 1)(s + 2) = s^2 + 3s + 2
  let pa = Polynomial::new([1.0, 1.0])
  let pb = Polynomial::new([1.0, 2.0])
  let prod = pa * pb
  assert_eq(prod.coeffs, [1.0, 3.0, 2.0])
}

test "polynomial derivative" {
  // (3s^3 + 2s^2 - 5s + 4)' = 9s^2 + 4s - 5
  let p = Polynomial::new([3.0, 2.0, -5.0, 4.0])
  let d = p.derivative()
  assert_eq(d.coeffs, [9.0, 4.0, -5.0])
}

test "polynomial roots degree 1 and 2" {
  // 2s + 6 = 0 -> s = -3
  let p1 = Polynomial::new([2.0, 6.0])
  let r1 = p1.roots()
  assert_eq(r1.length(), 1)
  assert_true(r1[0].approx_eq(Complex::from_real(-3.0)))

  // s^2 - 5s + 6 = 0 -> s = 2, 3
  let p2 = Polynomial::new([1.0, -5.0, 6.0])
  let r2 = p2.roots()
  assert_eq(r2.length(), 2)
  let found2 = r2[0].approx_eq(Complex::from_real(2.0)) || r2[1].approx_eq(Complex::from_real(2.0))
  let found3 = r2[0].approx_eq(Complex::from_real(3.0)) || r2[1].approx_eq(Complex::from_real(3.0))
  assert_true(found2 && found3)

  // s^2 + 2s + 5 = 0 -> s = -1 +/- 2i
  let p3 = Polynomial::new([1.0, 2.0, 5.0])
  let r3 = p3.roots()
  assert_eq(r3.length(), 2)
  assert_true(r3[0].approx_eq(Complex::new(-1.0, 2.0)) || r3[0].approx_eq(Complex::new(-1.0, -2.0)))
  assert_true(r3[1].approx_eq(Complex::new(-1.0, 2.0)) || r3[1].approx_eq(Complex::new(-1.0, -2.0)))
}

test "polynomial roots degree 3 and 4 Durand-Kerner" {
  // (s - 1)(s - 2)(s - 3) = s^3 - 6s^2 + 11s - 6
  let p_cube = Polynomial::new([1.0, -6.0, 11.0, -6.0])
  let roots3 = p_cube.roots()
  assert_eq(roots3.length(), 3)
  for r in roots3 {
    // Each root should evaluate p_cube(r) to ~ 0
    let ev = p_cube.eval_complex(r)
    assert_true(ev.abs() < 1.0e-7)
  }

  // (s^2 + 1)(s^2 + 4) = s^4 + 5s^2 + 4 -> roots: +/- i, +/- 2i
  let p4 = Polynomial::new([1.0, 0.0, 5.0, 0.0, 4.0])
  let roots4 = p4.roots()
  assert_eq(roots4.length(), 4)
  for r in roots4 {
    let ev = p4.eval_complex(r)
    assert_true(ev.abs() < 1.0e-7)
  }
}
