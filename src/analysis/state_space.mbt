///|
/// State-space structural properties: Kalman controllability, observability, and minimality.
pub(all) struct StateSpaceAnalysis {
  order : Int
  controllability_rank : Int
  observability_rank : Int
  is_controllable : Bool
  is_observable : Bool
  is_minimal : Bool
} derive(Eq)

pub extend StateSpaceAnalysis with Eq::{equal, not_equal}
pub extend StateSpaceAnalysis with Show::{output}

///|
/// Compute numerical rank of a matrix via Gaussian elimination with complete row/column pivoting.
pub fn matrix_rank(m : @num.Matrix, tol~ : Double = 1.0e-9) -> Int {
  let rows = m.rows
  let cols = m.cols
  let work = Array::make(rows * cols, 0.0)
  for i = 0; i < rows * cols; i = i + 1 {
    work[i] = m.data[i]
  }

  let mut rank = 0
  let mut r = 0
  let mut c = 0
  while r < rows && c < cols {
    // Find pivot in submatrix [r..rows, c]
    let mut max_r = r
    let mut max_val = work[r * cols + c].abs()
    for i = r + 1; i < rows; i = i + 1 {
      let v = work[i * cols + c].abs()
      if v > max_val {
        max_val = v
        max_r = i
      }
    }

    if max_val <= tol {
      c = c + 1
      continue
    }

    // Swap row r and max_r
    if max_r != r {
      for j = c; j < cols; j = j + 1 {
        let tmp = work[r * cols + j]
        work[r * cols + j] = work[max_r * cols + j]
        work[max_r * cols + j] = tmp
      }
    }

    // Eliminate below
    let pivot = work[r * cols + c]
    for i = r + 1; i < rows; i = i + 1 {
      let factor = work[i * cols + c] / pivot
      work[i * cols + c] = 0.0
      for j = c + 1; j < cols; j = j + 1 {
        work[i * cols + j] = work[i * cols + j] - factor * work[r * cols + j]
      }
    }

    rank = rank + 1
    r = r + 1
    c = c + 1
  }
  rank
}

///|
/// Construct the n x n Kalman controllability matrix C_mat = [B, A*B, A^2*B, ..., A^(n-1)*B]
/// for a single-input state-space pair (A: n x n, B: n x 1).
pub fn controllability_matrix(
  a : @num.Matrix,
  b : @num.Matrix
) -> Result[@num.Matrix, String] {
  let n = a.rows
  if a.cols != n {
    return Err("State matrix A must be square")
  }
  if b.rows != n || b.cols != 1 {
    return Err("Input matrix B must be n x 1")
  }

  let c_mat = @num.Matrix::zeros(n, n)
  let mut col_vec = b
  for col = 0; col < n; col = col + 1 {
    for row = 0; row < n; row = row + 1 {
      c_mat.set(row, col, col_vec.get(row, 0))
    }
    if col < n - 1 {
      col_vec = a * col_vec
    }
  }
  Ok(c_mat)
}

///|
/// Construct the n x n Kalman observability matrix O_mat = [C; C*A; C*A^2; ...; C*A^(n-1)]
/// for a single-output state-space pair (A: n x n, C: 1 x n).
pub fn observability_matrix(
  a : @num.Matrix,
  c : @num.Matrix
) -> Result[@num.Matrix, String] {
  let n = a.rows
  if a.cols != n {
    return Err("State matrix A must be square")
  }
  if c.rows != 1 || c.cols != n {
    return Err("Output matrix C must be 1 x n")
  }

  let o_mat = @num.Matrix::zeros(n, n)
  let mut row_vec = c
  for row = 0; row < n; row = row + 1 {
    for col = 0; col < n; col = col + 1 {
      o_mat.set(row, col, row_vec.get(0, col))
    }
    if row < n - 1 {
      row_vec = row_vec * a
    }
  }
  Ok(o_mat)
}

///|
/// Convert a proper SISO transfer function into explicit state-space matrices (A, B, C)
/// in controllable canonical form and evaluate Kalman controllability/observability ranks.
pub fn analyze_transfer_function_state_space(
  tf : @tf.TransferFunction
) -> Result[(@num.Matrix, @num.Matrix, @num.Matrix, StateSpaceAnalysis), String] {
  let n = tf.den.degree()
  if n <= 0 {
    return Err("System order must be at least 1")
  }
  if tf.num.degree() > n {
    return Err("Transfer function must be proper")
  }

  let a = @num.Matrix::zeros(n, n)
  for i = 0; i < n - 1; i = i + 1 {
    a.set(i, i + 1, 1.0)
  }
  for j = 0; j < n; j = j + 1 {
    let a_j = tf.den.coeffs[n - j]
    a.set(n - 1, j, -a_j)
  }

  let b = @num.Matrix::zeros(n, 1)
  b.set(n - 1, 0, 1.0)

  let c = @num.Matrix::zeros(1, n)
  let num_len = tf.num.coeffs.length()
  let d_gain = if tf.num.degree() == n { tf.num.coeffs[0] } else { 0.0 }
  for j = 0; j < n; j = j + 1 {
    let idx = num_len - 1 - j
    let b_j = if idx >= 0 { tf.num.coeffs[idx] } else { 0.0 }
    let a_j = tf.den.coeffs[n - j]
    c.set(0, j, b_j - d_gain * a_j)
  }

  let c_mat = match controllability_matrix(a, b) {
    Ok(m) => m
    Err(e) => return Err(e)
  }
  let o_mat = match observability_matrix(a, c) {
    Ok(m) => m
    Err(e) => return Err(e)
  }

  let c_rank = matrix_rank(c_mat)
  let o_rank = matrix_rank(o_mat)
  let is_ctrl = c_rank == n
  let is_obs = o_rank == n

  Ok((
    a,
    b,
    c,
    {
      order: n,
      controllability_rank: c_rank,
      observability_rank: o_rank,
      is_controllable: is_ctrl,
      is_observable: is_obs,
      is_minimal: is_ctrl && is_obs,
    },
  ))
}

///|
/// Ackermann's pole placement formula: K = [0, ..., 0, 1] * C_mat^(-1) * alpha_c(A)
/// Computes full-state feedback gain row vector K = [k_1, ..., k_n] such that the eigenvalues
/// of (A - B * K) match the roots of desired_char_poly: s^n + alpha_{n-1} s^(n-1) + ... + alpha_0.
pub fn ackermann_pole_placement(
  a : @num.Matrix,
  b : @num.Matrix,
  desired_char_poly : @num.Polynomial
) -> Result[Array[Double], String] {
  let n = a.rows
  if desired_char_poly.degree() != n {
    return Err("Desired characteristic polynomial degree must equal state dimension n")
  }
  let c_mat = match controllability_matrix(a, b) {
    Ok(m) => m
    Err(e) => return Err(e)
  }
  let c_inv = match c_mat.inv() {
    Ok(inv) => inv
    Err(_) => return Err("System is not completely state controllable (C_mat is singular)")
  }

  // Normalize desired_char_poly to be monic
  let lead = desired_char_poly.coeffs[0]
  let poly = desired_char_poly.scale(1.0 / lead)

  // Evaluate matrix polynomial alpha_c(A) = A^n + alpha_{n-1} A^(n-1) + ... + alpha_0 I
  // using Horner's method: P = I; P = P * A + c_i * I
  let mut phi_a = @num.Matrix::identity(n)
  for i = 1; i <= n; i = i + 1 {
    let coeff = poly.coeffs[i]
    phi_a = phi_a * a + @num.Matrix::identity(n).scale(coeff)
  }

  // e_n^T * C_mat^(-1) * phi_a is the last row (index n - 1) of (c_inv * phi_a)
  let prod = c_inv * phi_a
  let gains = Array::make(n, 0.0)
  for j = 0; j < n; j = j + 1 {
    gains[j] = prod.get(n - 1, j)
  }
  Ok(gains)
}

///|
pub fn StateSpaceAnalysis::to_string(self : StateSpaceAnalysis) -> String {
  "StateSpaceAnalysis(order: \{self.order}, ctrl_rank: \{self.controllability_rank}, obs_rank: \{self.observability_rank}, minimal: \{self.is_minimal})"
}

///|
pub impl Show for StateSpaceAnalysis with output(self, logger) {
  logger.write_string(self.to_string())
}
