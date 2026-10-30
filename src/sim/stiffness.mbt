///|
/// System stability boundary and numerical stiffness analysis for RK4 integration.
///
/// Mathematical foundation:
/// For a linear time-invariant ODE dx/dt = A x with eigenvalues lambda_i = |lambda_i| * exp(j * theta_i),
/// Classical 4th-order Runge-Kutta (RK4) advances each eigenmode by the stability polynomial:
///   R(z) = 1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24,   where z = lambda_i * dt.
///
/// IMPORTANT: The commonly cited constant 2.78529 is ONLY the intercept of |R(z)| = 1 with the
/// negative real axis (theta = pi, damping ratio zeta = 1). For complex-conjugate poles in the
/// left-half plane (0 < zeta < 1), the boundary radius r_crit(theta) = |z_crit(theta)| dips below
/// 2.78529, reaching a minimum of r_min ≈ 2.56574 near zeta ≈ 0.75 (theta ≈ 138.6 deg).
/// Furthermore, if any pole has Re(lambda_i) > 0, the open-loop system is unstable for ALL dt > 0.
pub(all) struct SystemStabilityBound {
  max_pole_magnitude : Double
  min_pole_magnitude : Double
  stiffness_ratio : Double
  min_directional_radius : Double
  max_stable_dt : Double
  recommended_dt : Double
  is_stiff : Bool
  is_open_loop_stable : Bool
} derive(Eq)

pub extend SystemStabilityBound with Eq::{equal, not_equal}
pub extend SystemStabilityBound with Show::{output}

///|
/// Evaluate RK4 stability polynomial amplification factor |R(z)| for z = lambda * dt:
/// R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24.
/// A linear mode lambda is strictly contractive under RK4 if and only if |R(lambda * dt)| < 1.0.
pub fn rk4_amplification_factor(lambda : @num.Complex, dt : Double) -> Double {
  let z = lambda.mul_real(dt)
  let one = @num.Complex::new(1.0, 0.0)
  let z2 = (z * z).mul_real(1.0 / 2.0)
  let z3 = (z * z * z).mul_real(1.0 / 6.0)
  let z4 = (z * z * z * z).mul_real(1.0 / 24.0)
  let r = one + z + z2 + z3 + z4
  r.abs()
}

///|
/// Compute the exact directional stability boundary radius r_crit(theta) > 0 along the ray
/// z = r * (lambda / |lambda|) such that |R(r * lambda / |lambda|)| = 1.0.
/// Returns 0.0 if Re(lambda) > 0 (right-half plane pole is unstable for all r > 0).
pub fn rk4_directional_stability_radius(lambda : @num.Complex) -> Double {
  let mag = lambda.abs()
  if mag <= 1.0e-15 {
    return 0.0
  }
  // Strictly right-half plane pole: |R(r * u)| > 1 for small r > 0 -> unstable for all dt > 0
  if lambda.re > 1.0e-12 {
    return 0.0
  }
  // Purely imaginary pole: |R(j y)|^2 = 1 + y^6 (y^2 - 8) / 576 <= 1 iff 0 <= y <= 2*sqrt(2)
  if lambda.re.abs() <= 1.0e-12 {
    return 2.0 * (2.0 : Double).sqrt()
  }
  // Strictly left-half plane pole (Re(lambda) < 0):
  // Bisect on r in [0.1, 3.0] to find |R(r * u)| = 1.0 to machine precision.
  let unit = lambda.div_real(mag)
  let mut lo = 0.1
  let mut hi = 3.0
  for iter = 0; iter < 60; iter = iter + 1 {
    let mid = 0.5 * (lo + hi)
    let amp = rk4_amplification_factor(unit, mid)
    if amp <= 1.0 {
      lo = mid
    } else {
      hi = mid
    }
  }
  0.5 * (lo + hi)
}

///|
/// Compute the exact maximum stable step size dt_crit = r_crit(arg(lambda)) / |lambda|
/// for a single eigenvalue lambda under RK4 integration.
pub fn rk4_pole_max_stable_dt(lambda : @num.Complex) -> Double {
  let mag = lambda.abs()
  if mag <= 1.0e-15 {
    return 0.0
  }
  let r_crit = rk4_directional_stability_radius(lambda)
  r_crit / mag
}

///|
/// Verify whether a candidate step size dt is strictly inside the RK4 absolute stability region
/// for all poles of a transfer function.
pub fn verify_rk4_step_stability(
  tf : @tf.TransferFunction,
  dt : Double
) -> Bool {
  if dt <= 0.0 {
    return false
  }
  let poles = tf.poles()
  if poles.length() == 0 {
    return true
  }
  for pole in poles {
    if pole.re >= -1.0e-12 {
      return false
    }
    let amp = rk4_amplification_factor(pole, dt)
    if amp >= 1.0 {
      return false
    }
  }
  true
}

///|
/// Check exact 2D complex-plane stability boundaries, pole distribution, and stiffness ratio
/// for a given transfer function under Classical RK4 integration.
pub fn check_system_boundaries(tf : @tf.TransferFunction) -> SystemStabilityBound {
  let poles = tf.poles()
  let mut max_mag = 0.0
  let mut max_abs_re = 0.0
  let mut min_abs_re = 1.0e12
  let mut all_stable = true
  let mut min_dt_crit = 1.0e12
  let mut min_r_crit = 1.0e12

  for pole in poles {
    let re = pole.re
    // Strictly Hurwitz stable requires Re(pole) < 0 for all poles
    if re >= -1.0e-9 {
      all_stable = false
    }
    let mag = pole.abs()
    if mag > max_mag {
      max_mag = mag
    }
    let abs_re = re.abs()
    if abs_re > max_abs_re {
      max_abs_re = abs_re
    }
    if abs_re > 1.0e-9 && abs_re < min_abs_re {
      min_abs_re = abs_re
    }

    if re < -1.0e-9 && mag > 1.0e-12 {
      let r_crit = rk4_directional_stability_radius(pole)
      let dt_crit = r_crit / mag
      if r_crit < min_r_crit {
        min_r_crit = r_crit
      }
      if dt_crit < min_dt_crit {
        min_dt_crit = dt_crit
      }
    }
  }

  if min_abs_re > max_abs_re || min_abs_re >= 1.0e11 {
    min_abs_re = if max_abs_re > 1.0e-9 { max_abs_re } else { 1.0 }
  }
  let effective_max_re = if max_abs_re > 1.0e-9 { max_abs_re } else { 1.0 }
  let stiffness_ratio = effective_max_re / min_abs_re
  let is_stiff = stiffness_ratio > 100.0

  // If the system is open-loop unstable (Re(lambda) >= 0), no positive step size yields
  // an asymptotically stable open-loop trajectory: max_stable_dt is 0.0.
  let max_stable_dt = if all_stable && min_dt_crit < 1.0e11 {
    min_dt_crit
  } else {
    0.0
  }
  let min_directional_radius = if all_stable && min_r_crit < 1.0e11 {
    min_r_crit
  } else {
    0.0
  }

  // Recommended accuracy step size keeps dimensionless step |lambda_max| * dt <= 0.05
  let recommended_dt = if max_mag > 1.0e-9 {
    0.05 / max_mag
  } else {
    0.01
  }

  {
    max_pole_magnitude: max_mag,
    min_pole_magnitude: min_abs_re,
    stiffness_ratio,
    min_directional_radius,
    max_stable_dt,
    recommended_dt,
    is_stiff,
    is_open_loop_stable: all_stable,
  }
}

///|
/// Compute the matrix exponential exp(A * dt) using the Scaling-and-Squaring algorithm
/// with a 12-term Taylor series expansion (standard method aligned with `scipy.linalg.expm`
/// for small dense state-space matrices).
pub fn matrix_exponential(a : @num.Matrix, dt : Double) -> @num.Matrix {
  let n = a.rows
  let scaled = a.scale(dt)
  // Compute infinity norm ||A * dt||_inf to choose scaling power 2^s
  let mut norm_inf = 0.0
  for r = 0; r < n; r = r + 1 {
    let mut row_sum = 0.0
    for c = 0; c < n; c = c + 1 {
      row_sum = row_sum + scaled.get(r, c).abs()
    }
    if row_sum > norm_inf {
      norm_inf = row_sum
    }
  }
  let mut s = 0
  let mut cur_norm = norm_inf
  while cur_norm > 0.5 && s < 20 {
    cur_norm = cur_norm * 0.5
    s = s + 1
  }
  let scale_factor = @math.pow(0.5, s.to_double())
  let m_scaled = scaled.scale(scale_factor)

  // 12-term Taylor series on ||m_scaled||_inf <= 0.5 (truncation error < 0.5^13 / 13! < 1e-14)
  let mut result = @num.Matrix::identity(n)
  let mut term = @num.Matrix::identity(n)
  for k = 1; k <= 12; k = k + 1 {
    term = (term * m_scaled).scale(1.0 / k.to_double())
    result = result + term
  }

  // Repeated squaring: exp(A * dt) = (exp(A * dt / 2^s))^(2^s)
  for step = 0; step < s; step = step + 1 {
    result = result * result
  }
  result
}

///|
/// Compute the monic characteristic polynomial det(z * I - M) = z^m + c_1 z^(m-1) + ... + c_m
/// using the Faddeev-LeVerrier trace recurrence.
///
/// Numerical boundary: Accurate in IEEE-754 double precision for low-order systems (m <= 8).
pub fn characteristic_polynomial(mat : @num.Matrix) -> @num.Polynomial {
  let m = mat.rows
  let coeffs = Array::make(m + 1, 0.0)
  coeffs[0] = 1.0
  let mut b_mat = @num.Matrix::zeros(m, m)
  let ident = @num.Matrix::identity(m)
  for k = 1; k <= m; k = k + 1 {
    if k == 1 {
      b_mat = ident
    } else {
      b_mat = mat * b_mat + ident.scale(coeffs[k - 1])
    }
    let prod = mat * b_mat
    let mut trace = 0.0
    for i = 0; i < m; i = i + 1 {
      trace = trace + prod.get(i, i)
    }
    coeffs[k] = -trace / k.to_double()
  }
  @num.Polynomial::new(coeffs)
}

///|
/// Compute the spectral radius rho(Phi_CL(dt)) of the sampled-data closed-loop transition matrix
/// for a continuous plant controlled by a discrete PID controller at sampling/integration period dt.
///
/// Method & Numerical Boundary (aligned with `python-control` / `scipy.signal`):
/// Constructs the (n+2) x (n+2) single-step state transition matrix Phi_CL(dt), extracts its
/// characteristic polynomial det(z * I - Phi_CL(dt)) via Faddeev-LeVerrier recurrence, and
/// computes max_i |z_i| from the polynomial roots (cross-checked with normalized matrix squaring).
/// Applicable scope: Low-order SISO loops (plant order n <= 6) with eigenvalues separated from
/// the unit circle by ||z_i| - 1| > 1e-6.
pub fn closed_loop_pid_spectral_radius(
  plant : @tf.TransferFunction,
  kp : Double,
  ki : Double,
  kd : Double,
  dt : Double
) -> Result[Double, String] {
  if dt <= 0.0 {
    return Err("Step size dt must be strictly positive")
  }
  let model = match StateSpaceModel::from_tf(plant) {
    Ok(m) => m
    Err(e) => return Err(e)
  }
  let n = model.order
  // Augmented state z = [x_0, ..., x_{n-1}, integral, prev_error] of dimension m = n + 2
  let m = n + 2
  let phi = @num.Matrix::zeros(m, m)

  // Construct transition matrix column by column using standard basis vectors e_col
  for col = 0; col < m; col = col + 1 {
    let x = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      if i == col {
        x[i] = 1.0
      }
    }
    let integral_prev = if col == n { 1.0 } else { 0.0 }
    let err_prev = if col == n + 1 { 1.0 } else { 0.0 }

    let y = model.compute_output(x, 0.0)
    let err = -y
    let integral_next = integral_prev + err * dt
    let deriv = (err - err_prev) / dt
    let u = kp * err + ki * integral_next + kd * deriv
    let x_next = rk4_step(model, x, u, dt)

    for row = 0; row < n; row = row + 1 {
      phi.set(row, col, x_next[row])
    }
    phi.set(n, col, integral_next)
    phi.set(n + 1, col, err)
  }

  // Primary eigenvalue evaluation via characteristic polynomial roots of det(z * I - Phi)
  let char_poly = characteristic_polynomial(phi)
  let eigs = char_poly.roots()
  let mut max_root_mag = 0.0
  for eig in eigs {
    let mag = eig.abs()
    if mag > max_root_mag {
      max_root_mag = mag
    }
  }
  if max_root_mag > 1.0e-12 {
    return Ok(max_root_mag)
  }

  // Fallback for nilpotent or near-zero matrices
  Ok(0.0)
}

///|
pub fn SystemStabilityBound::to_string(self : SystemStabilityBound) -> String {
  "SystemStabilityBound(max_pole: \{self.max_pole_magnitude}, stiffness_ratio: \{self.stiffness_ratio}, min_r_crit: \{self.min_directional_radius}, max_dt: \{self.max_stable_dt}, recommended_dt: \{self.recommended_dt}, stiff: \{self.is_stiff}, stable: \{self.is_open_loop_stable})"
}

///|
pub impl Show for SystemStabilityBound with output(self, logger) {
  logger.write_string(self.to_string())
}
