///|
/// Dense 2D matrix stored in row-major order: `data[r * cols + c]`
pub(all) struct Matrix {
  rows : Int
  cols : Int
  data : Array[Double]
} derive(Eq)

pub extend Matrix with Eq::{equal, not_equal}
pub extend Matrix with Show::{output}
pub extend Matrix with Add::{add}
pub extend Matrix with Sub::{sub}
pub extend Matrix with Mul::{mul}

///|
pub fn Matrix::new(rows : Int, cols : Int, data : Array[Double]) -> Matrix {
  { rows, cols, data }
}

///|
pub fn Matrix::zeros(rows : Int, cols : Int) -> Matrix {
  { rows, cols, data: Array::make(rows * cols, 0.0) }
}

///|
pub fn Matrix::identity(n : Int) -> Matrix {
  let data = Array::make(n * n, 0.0)
  for i = 0; i < n; i = i + 1 {
    data[i * n + i] = 1.0
  }
  { rows: n, cols: n, data }
}

///|
pub fn Matrix::get(self : Matrix, r : Int, c : Int) -> Double {
  self.data[r * self.cols + c]
}

///|
pub fn Matrix::set(self : Matrix, r : Int, c : Int, val : Double) -> Unit {
  self.data[r * self.cols + c] = val
}

///|
pub impl Add for Matrix with add(self, other) {
  let len = self.rows * self.cols
  let res = Array::make(len, 0.0)
  for i = 0; i < len; i = i + 1 {
    res[i] = self.data[i] + other.data[i]
  }
  { rows: self.rows, cols: self.cols, data: res }
}

///|
pub impl Sub for Matrix with sub(self, other) {
  let len = self.rows * self.cols
  let res = Array::make(len, 0.0)
  for i = 0; i < len; i = i + 1 {
    res[i] = self.data[i] - other.data[i]
  }
  { rows: self.rows, cols: self.cols, data: res }
}

///|
pub impl Mul for Matrix with mul(self, other) {
  let r1 = self.rows
  let c1 = self.cols
  let c2 = other.cols
  let res = Array::make(r1 * c2, 0.0)
  for i = 0; i < r1; i = i + 1 {
    for k = 0; k < c1; k = k + 1 {
      let a_ik = self.data[i * c1 + k]
      for j = 0; j < c2; j = j + 1 {
        let idx = i * c2 + j
        res[idx] = res[idx] + a_ik * other.data[k * c2 + j]
      }
    }
  }
  { rows: r1, cols: c2, data: res }
}

///|
pub fn Matrix::scale(self : Matrix, factor : Double) -> Matrix {
  let len = self.rows * self.cols
  let res = Array::make(len, 0.0)
  for i = 0; i < len; i = i + 1 {
    res[i] = self.data[i] * factor
  }
  { rows: self.rows, cols: self.cols, data: res }
}

///|
pub fn Matrix::transpose(self : Matrix) -> Matrix {
  let res = Array::make(self.rows * self.cols, 0.0)
  for r = 0; r < self.rows; r = r + 1 {
    for c = 0; c < self.cols; c = c + 1 {
      res[c * self.rows + r] = self.data[r * self.cols + c]
    }
  }
  { rows: self.cols, cols: self.rows, data: res }
}

///|
/// Inverse using Gauss-Jordan elimination with partial pivoting.
pub fn Matrix::inv(self : Matrix) -> Result[Matrix, String] {
  if self.rows != self.cols {
    return Err("Matrix must be square to invert")
  }
  let n = self.rows
  // Augmented matrix [A | I]
  let aug = Array::make(n * (2 * n), 0.0)
  let stride = 2 * n
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      aug[i * stride + j] = self.data[i * n + j]
    }
    aug[i * stride + n + i] = 1.0
  }

  for col = 0; col < n; col = col + 1 {
    // Find pivot
    let mut max_row = col
    let mut max_val = aug[col * stride + col].abs()
    for row = col + 1; row < n; row = row + 1 {
      let val = aug[row * stride + col].abs()
      if val > max_val {
        max_val = val
        max_row = row
      }
    }

    if max_val < 1.0e-12 {
      return Err("Matrix is singular (non-invertible)")
    }

    // Swap pivot row
    if max_row != col {
      for j = 0; j < stride; j = j + 1 {
        let tmp = aug[col * stride + j]
        aug[col * stride + j] = aug[max_row * stride + j]
        aug[max_row * stride + j] = tmp
      }
    }

    // Normalize pivot row
    let pivot = aug[col * stride + col]
    for j = 0; j < stride; j = j + 1 {
      aug[col * stride + j] = aug[col * stride + j] / pivot
    }

    // Eliminate other rows
    for row = 0; row < n; row = row + 1 {
      if row != col {
        let factor = aug[row * stride + col]
        for j = 0; j < stride; j = j + 1 {
          aug[row * stride + j] = aug[row * stride + j] - factor * aug[col * stride + j]
        }
      }
    }
  }

  // Extract right half
  let inv_data = Array::make(n * n, 0.0)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      inv_data[i * n + j] = aug[i * stride + n + j]
    }
  }
  Ok({ rows: n, cols: n, data: inv_data })
}

///|
pub fn Matrix::to_string(self : Matrix) -> String {
  let buf = StringBuilder()
  buf.write_string("[\n")
  for r = 0; r < self.rows; r = r + 1 {
    buf.write_string("  [")
    for c = 0; c < self.cols; c = c + 1 {
      if c > 0 {
        buf.write_string(", ")
      }
      buf.write_string("\{self.data[r * self.cols + c]}")
    }
    buf.write_string("]\n")
  }
  buf.write_string("]")
  buf.to_string()
}

///|
pub impl Show for Matrix with output(self, logger) {
  logger.write_string(self.to_string())
}
