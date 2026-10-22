test "transfer function creation and normalization" {
  let g = TransferFunction::from_coeffs([4.0], [2.0, 6.0]).unwrap()
  // 4 / (2s + 6) normalized to monic den -> 2 / (s + 3)
  assert_eq(g.num.coeffs, [2.0])
  assert_eq(g.den.coeffs, [1.0, 3.0])
  assert_eq(g.order(), 1)
}

test "transfer function series and parallel" {
  let g1 = TransferFunction::from_coeffs([1.0], [1.0, 1.0]).unwrap() // 1 / (s + 1)
  let g2 = TransferFunction::from_coeffs([2.0], [1.0, 2.0]).unwrap() // 2 / (s + 2)

  // Series: 2 / (s^2 + 3s + 2)
  let g_series = g1 * g2
  assert_eq(g_series.num.coeffs, [2.0])
  assert_eq(g_series.den.coeffs, [1.0, 3.0, 2.0])

  // Parallel: (3s + 4) / (s^2 + 3s + 2)
  let g_parallel = g1 + g2
  assert_eq(g_parallel.num.coeffs, [3.0, 4.0])
  assert_eq(g_parallel.den.coeffs, [1.0, 3.0, 2.0])
}

test "transfer function feedback connection" {
  // G = 10 / (s^2 + 2s), unity negative feedback -> G_cl = 10 / (s^2 + 2s + 10)
  let g = TransferFunction::from_coeffs([10.0], [1.0, 2.0, 0.0]).unwrap()
  let g_cl = g.unity_feedback()

  assert_eq(g_cl.num.coeffs, [10.0])
  assert_eq(g_cl.den.coeffs, [1.0, 2.0, 10.0])
}

test "transfer function poles zeros and dc gain" {
  // G = (2s + 6) / (s^2 + 2s + 5)
  // Zero: s = -3
  // Poles: s = -1 +/- 2i
  let g = TransferFunction::from_coeffs([2.0, 6.0], [1.0, 2.0, 5.0]).unwrap()

  let zs = g.zeros()
  assert_eq(zs.length(), 1)
  assert_true(zs[0].approx_eq(@num.Complex::from_real(-3.0)))

  let ps = g.poles()
  assert_eq(ps.length(), 2)
  assert_true(ps[0].approx_eq(@num.Complex::new(-1.0, 2.0)) || ps[0].approx_eq(@num.Complex::new(-1.0, -2.0)))
  assert_true(ps[1].approx_eq(@num.Complex::new(-1.0, 2.0)) || ps[1].approx_eq(@num.Complex::new(-1.0, -2.0)))

  // DC gain: 6 / 5 = 1.2
  match g.dc_gain() {
    Some(k) => assert_true((k - 1.2).abs() < 1.0e-9)
    None => assert_true(false)
  }

  // Integrator DC gain should be None
  let integ = TransferFunction::from_coeffs([1.0], [1.0, 0.0]).unwrap()
  assert_true(integ.dc_gain() is None)
}
