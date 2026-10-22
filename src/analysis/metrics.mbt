///|
/// Performance metrics for underdamped second-order step response
pub(all) struct PerformanceMetrics {
  overshoot_fraction : Double
  overshoot_percent : Double
  peak_time : Double
  rise_time : Double
  settling_time_2pct : Double
  settling_time_5pct : Double
} derive(Eq)

pub extend PerformanceMetrics with Eq::{equal, not_equal}
pub extend PerformanceMetrics with Show::{output}

///|
/// Calculate standard second-order performance metrics.
pub fn calculate_second_order_metrics(
  omega_n : Double,
  zeta : Double
) -> Result[PerformanceMetrics, String] {
  if omega_n <= 0.0 {
    return Err("omega_n must be strictly positive")
  }
  if zeta <= 0.0 || zeta >= 1.0 {
    return Err("zeta must be in (0, 1) for underdamped performance metrics")
  }
  let sqrt_1_minus_zeta2 = (1.0 - zeta * zeta).sqrt()
  let omega_d = omega_n * sqrt_1_minus_zeta2
  let phi = @math.acos(zeta)

  // Mp = exp(-pi * zeta / sqrt(1 - zeta^2))
  let mp = @math.exp(-@math.PI * zeta / sqrt_1_minus_zeta2)
  let mp_pct = mp * 100.0

  // tp = pi / omega_d
  let tp = @math.PI / omega_d

  // tr = (pi - phi) / omega_d
  let tr = (@math.PI - phi) / omega_d

  // ts_2pct ≈ 4 / (zeta * omega_n)
  let ts_2pct = 4.0 / (zeta * omega_n)

  // ts_5pct ≈ 3 / (zeta * omega_n)
  let ts_5pct = 3.0 / (zeta * omega_n)

  Ok({
    overshoot_fraction: mp,
    overshoot_percent: mp_pct,
    peak_time: tp,
    rise_time: tr,
    settling_time_2pct: ts_2pct,
    settling_time_5pct: ts_5pct,
  })
}

///|
pub fn PerformanceMetrics::to_string(self : PerformanceMetrics) -> String {
  "PerformanceMetrics { Mp: \{self.overshoot_percent}%, tp: \{self.peak_time}s, tr: \{self.rise_time}s, ts(2%): \{self.settling_time_2pct}s }"
}

///|
pub impl Show for PerformanceMetrics with output(self, logger) {
  logger.write_string(self.to_string())
}

///|
pub(all) enum InputType {
  Step
  Ramp
  Parabolic
} derive(Eq)

pub extend InputType with Eq::{equal, not_equal}

///|
/// Calculate steady-state error e_ss for unity feedback system with open-loop transfer function G(s).
/// Uses Final Value Theorem based on system type (number of pure integrators at origin).
pub fn steady_state_error(
  open_loop : @tf.TransferFunction,
  input : InputType
) -> Result[Double, String] {
  let den = open_loop.den.coeffs
  let num = open_loop.num.coeffs
  let den_len = den.length()
  let num_len = num.length()

  // Count trailing zeros in denominator to determine system type N
  let mut sys_type = 0
  let mut idx = den_len - 1
  while idx >= 0 && den[idx].abs() < 1.0e-12 {
    sys_type = sys_type + 1
    idx = idx - 1
  }

  // Value of non-zero trailing coefficient
  if idx < 0 {
    return Err("Invalid denominator: all coefficients are zero")
  }
  let den_constant = den[idx]
  let num_constant = num[num_len - 1]
  let gain = num_constant / den_constant

  match input {
    Step =>
      if sys_type == 0 {
        // e_ss = 1 / (1 + Kp), Kp = gain
        Ok(1.0 / (1.0 + gain))
      } else {
        // Type >= 1 has zero steady state error to step
        Ok(0.0)
      }
    Ramp =>
      if sys_type == 0 {
        // Infinite error for Type 0
        Ok(@double.infinity)
      } else if sys_type == 1 {
        // e_ss = 1 / Kv, Kv = gain
        Ok(1.0 / gain)
      } else {
        // Type >= 2 has zero steady state error to ramp
        Ok(0.0)
      }
    Parabolic =>
      if sys_type < 2 {
        Ok(@double.infinity)
      } else if sys_type == 2 {
        // e_ss = 1 / Ka, Ka = gain
        Ok(1.0 / gain)
      } else {
        Ok(0.0)
      }
  }
}
