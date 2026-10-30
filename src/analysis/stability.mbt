///|
/// Routh Array representation for stability analysis
pub(all) struct RouthArray {
  rows : Array[Array[Double]]
  sign_changes : Int
  is_stable : Bool
} derive(Eq)

pub extend RouthArray with Eq::{equal, not_equal}
pub extend RouthArray with Show::{output}

///|
/// Construct Routh-Hurwitz array and evaluate system stability from characteristic polynomial.
pub fn routh_hurwitz(char_poly : @num.Polynomial) -> Result[RouthArray, String] {
  let deg = char_poly.degree()
  if deg <= 0 {
    return Err("Characteristic polynomial degree must be at least 1")
  }

  // Normalize so leading coefficient is positive
  let poly = if char_poly.coeffs[0] < 0.0 {
    char_poly.scale(-1.0)
  } else {
    char_poly
  }

  let num_rows = deg + 1
  let num_cols = deg / 2 + 1
  let table = Array::make(num_rows, Array::make(num_cols, 0.0))
  for i = 0; i < num_rows; i = i + 1 {
    table[i] = Array::make(num_cols, 0.0)
  }

  // Fill Row 0 (s^n)
  let mut c = 0
  let mut p = 0
  while p <= deg {
    table[0][c] = poly.coeffs[p]
    c = c + 1
    p = p + 2
  }

  // Fill Row 1 (s^(n-1))
  c = 0
  p = 1
  while p <= deg {
    table[1][c] = poly.coeffs[p]
    c = c + 1
    p = p + 2
  }

  // Necessary condition for strict Hurwitz stability: all polynomial coefficients must be > 0
  let mut all_coeffs_positive = true
  for c_val in poly.coeffs {
    if c_val <= 1.0e-12 {
      all_coeffs_positive = false
    }
  }

  let mut had_all_zero_row = false
  let mut had_zero_pivot = false

  // Compute remaining rows
  for i = 2; i < num_rows; i = i + 1 {
    // Check if previous row first element is near zero
    let prev_first = table[i - 1][0]
    if prev_first.abs() < 1.0e-9 {
      had_zero_pivot = true
      // Check if entire previous row is zero
      let mut all_zero = true
      for j = 0; j < num_cols; j = j + 1 {
        if table[i - 1][j].abs() >= 1.0e-9 {
          all_zero = false
          break
        }
      }

      if all_zero {
        had_all_zero_row = true
        // Auxiliary polynomial from row i-2: power is k = deg - (i - 2)
        let k = deg - (i - 2)
        for j = 0; j < num_cols; j = j + 1 {
          let power = k - 2 * j
          if power > 0 {
            table[i - 1][j] = table[i - 2][j] * power.to_double()
          } else {
            table[i - 1][j] = 0.0
          }
        }
      } else {
        // Perturb first element with small epsilon
        table[i - 1][0] = 1.0e-6
      }
    }

    let pivot = table[i - 1][0]
    for j = 0; j < num_cols - 1; j = j + 1 {
      let a1 = table[i - 1][0]
      let b1 = table[i - 2][j + 1]
      let a2 = table[i - 2][0]
      let b2 = table[i - 1][j + 1]
      table[i][j] = (a1 * b1 - a2 * b2) / pivot
    }
  }

  if table[num_rows - 1][0].abs() < 1.0e-9 {
    had_zero_pivot = true
  }

  // Count sign changes in first column
  let mut sign_changes = 0
  for i = 0; i < num_rows - 1; i = i + 1 {
    let curr = table[i][0]
    let next = table[i + 1][0]
    if curr * next < 0.0 {
      sign_changes = sign_changes + 1
    }
  }

  let stable = sign_changes == 0 && !had_all_zero_row && !had_zero_pivot && all_coeffs_positive
  Ok({ rows: table, sign_changes, is_stable: stable })
}

///|
/// Check stability by computing all poles of the transfer function and verifying Re(p) < 0.
pub fn is_stable_by_poles(
  tf : @tf.TransferFunction,
  tolerance~ : Double = 1.0e-9
) -> Bool {
  let ps = tf.poles()
  for p in ps {
    if p.re >= -tolerance {
      return false
    }
  }
  true
}

///|
pub fn RouthArray::to_string(self : RouthArray) -> String {
  let buf = StringBuilder()
  buf.write_string("RouthArray {\n")
  for r in self.rows {
    buf.write_string("  ")
    for v in r {
      buf.write_string("\{v}\t")
    }
    buf.write_string("\n")
  }
  buf.write_string(
    "  sign_changes: \{self.sign_changes}, is_stable: \{self.is_stable}\n}",
  )
  buf.to_string()
}

///|
pub impl Show for RouthArray with output(self, logger) {
  logger.write_string(self.to_string())
}
