///|
/// System stability boundary and numerical stiffness analysis for RK4 integration.
pub(all) struct SystemStabilityBound {
  max_pole_magnitude : Double
  min_pole_magnitude : Double
  stiffness_ratio : Double
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
/// The numerical integration is stable if and only if |R(z)| <= 1.0.
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
/// Check stability boundaries, pole distribution, and calculate maximum stable step size dt
/// for a given transfer function under Classical RK4 integration.
pub fn check_system_boundaries(tf : @tf.TransferFunction) -> SystemStabilityBound {
  let poles = tf.poles()
  let mut max_mag = 0.0
  let mut min_mag = 1.0e12
  let mut max_re = -1.0e12
  let mut all_stable = true

  for pole in poles {
    let re = pole.re
    if re > max_re {
      max_re = re
    }
    if re > 1.0e-9 {
      all_stable = false
    }
    let mag = pole.abs()
    if mag > max_mag {
      max_mag = mag
    }
    let abs_re = re.abs()
    if abs_re > 1.0e-6 && abs_re < min_mag {
      min_mag = abs_re
    }
  }

  if min_mag > max_mag {
    min_mag = max_mag
  }
  if min_mag < 1.0e-6 {
    min_mag = 1.0e-6
  }
  if max_mag < 1.0e-6 {
    max_mag = 1.0
    min_mag = 1.0
  }

  let stiffness_ratio = max_mag / min_mag
  let is_stiff = stiffness_ratio > 100.0

  // RK4 absolute stability limit along negative real axis is ~2.7853
  let max_stable_dt = 2.785 / max_mag
  let recommended_dt = 0.25 / max_mag

  {
    max_pole_magnitude: max_mag,
    min_pole_magnitude: min_mag,
    stiffness_ratio,
    max_stable_dt,
    recommended_dt,
    is_stiff,
    is_open_loop_stable: all_stable,
  }
}

///|
pub fn SystemStabilityBound::to_string(self : SystemStabilityBound) -> String {
  "SystemStabilityBound(max_pole: \{self.max_pole_magnitude}, stiffness_ratio: \{self.stiffness_ratio}, max_dt: \{self.max_stable_dt}, recommended_dt: \{self.recommended_dt}, stiff: \{self.is_stiff}, stable: \{self.is_open_loop_stable})"
}

///|
pub impl Show for SystemStabilityBound with output(self, logger) {
  logger.write_string(self.to_string())
}
