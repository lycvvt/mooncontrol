test "matrix operations" {
  let m1 = Matrix::new(2, 2, [1.0, 2.0, 3.0, 4.0])
  let m2 = Matrix::new(2, 2, [2.0, 0.0, 1.0, 2.0])

  let sum = m1 + m2
  assert_eq(sum.data, [3.0, 2.0, 4.0, 6.0])

  // [1 2; 3 4] * [2 0; 1 2] = [4 4; 10 8]
  let prod = m1 * m2
  assert_eq(prod.data, [4.0, 4.0, 10.0, 8.0])

  let trans = m1.transpose()
  assert_eq(trans.data, [1.0, 3.0, 2.0, 4.0])
}

test "matrix inverse" {
  // [4 7; 2 6] -> det = 24 - 14 = 10 -> inv = 0.1 * [6 -7; -2 4] = [0.6 -0.7; -0.2 0.4]
  let m = Matrix::new(2, 2, [4.0, 7.0, 2.0, 6.0])
  match m.inv() {
    Ok(inv_m) => {
      let ident = m * inv_m
      assert_true((ident.get(0, 0) - 1.0).abs() < 1.0e-9)
      assert_true(ident.get(0, 1).abs() < 1.0e-9)
      assert_true(ident.get(1, 0).abs() < 1.0e-9)
      assert_true((ident.get(1, 1) - 1.0).abs() < 1.0e-9)
    }
    Err(_) => assert_true(false)
  }
}
