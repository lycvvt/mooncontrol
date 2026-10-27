///|
/// Quantitative frequency-domain stability and robustness margins for an open-loop LTI system L(s):
/// - gain_margin_db: Gain Margin (dB) at phase crossover frequency omega_pc (where phase = -180 deg)
/// - phase_margin_deg: Phase Margin (deg) at gain crossover frequency omega_gc (where |L(j*omega)| = 1)
/// - gain_crossover_freq: omega_gc (rad/s)
/// - phase_crossover_freq: omega_pc (rad/s)
/// - sensitivity_peak: M_s = max_omega |1 / (1 + L(j*omega))| (robustness indicator; M_s < 2.0 is desirable)
///
/// Note: Unlike CAIMEOX/symbit/symphysics/control which only samples raw Bode arrays for plotting,
/// this module numerically solves for crossover frequencies via bisection and extracts quantitative margins.
pub(all) struct StabilityMargins {
  gain_margin_db : Double
  phase_margin_deg : Double
  gain_crossover_freq : Double
  phase_crossover_freq : Double
  sensitivity_peak : Double
} derive(Eq)

pub extend StabilityMargins with Eq::{equal, not_equal}
pub extend StabilityMargins with Show::{output}

///|
/// Compute quantitative Gain Margin, Phase Margin, Crossover Frequencies, and Sensitivity Peak
/// for an open-loop transfer function L(s).
pub fn compute_stability_margins(
  open_loop : @tf.TransferFunction,
  omega_min~ : Double = 0.001,
  omega_max~ : Double = 1000.0,
  num_points~ : Int = 1200
) -> StabilityMargins {
  let n = if num_points > 50 { num_points } else { 50 }
  let log_min = @math.ln(omega_min)
  let log_max = @math.ln(omega_max)
  let d_log = (log_max - log_min) / (n - 1).to_double()

  let mut prev_w = omega_min
  let mut prev_val = open_loop.eval_complex(@num.Complex::new(0.0, prev_w))
  let mut prev_mag = prev_val.abs()

  // Initial unwrapped phase in radians at omega_min
  let mut unwrapped_phase = prev_val.arg()
  if unwrapped_phase > 0.01 {
    unwrapped_phase = unwrapped_phase - 2.0 * @math.PI
  }
  let mut prev_phase = unwrapped_phase

  let mut w_gc = 0.0
  let mut pm_deg = @double.infinity
  let mut w_pc = 0.0
  let mut gm_db = @double.infinity
  let mut max_ms = 0.0

  let one = @num.Complex::new(1.0, 0.0)

  for i = 1; i < n; i = i + 1 {
    let w = @math.exp(log_min + d_log * i.to_double())
    let val = open_loop.eval_complex(@num.Complex::new(0.0, w))
    let mag = val.abs()

    // Continuous phase increment via arg(val * conj(prev_val))
    let conj_prev = @num.Complex::new(prev_val.re, -prev_val.im)
    let rel = val * conj_prev
    let d_phi = rel.arg()
    unwrapped_phase = prev_phase + d_phi

    // Sensitivity |1 / (1 + L(jw))|
    let denom_s = (one + val).abs()
    if denom_s > 1.0e-12 {
      let ms = 1.0 / denom_s
      if ms > max_ms {
        max_ms = ms
      }
    }

    // Check Gain Crossover: |L(jw)| crosses 1.0
    if (prev_mag >= 1.0 && mag < 1.0) || (prev_mag <= 1.0 && mag > 1.0) {
      let mut lo = prev_w
      let mut hi = w
      let mut lo_mag = prev_mag
      for iter = 0; iter < 50; iter = iter + 1 {
        let mid = 0.5 * (lo + hi)
        let mid_mag = open_loop.eval_complex(@num.Complex::new(0.0, mid)).abs()
        if (lo_mag >= 1.0 && mid_mag >= 1.0) || (lo_mag <= 1.0 && mid_mag <= 1.0) {
          lo = mid
          lo_mag = mid_mag
        } else {
          hi = mid
        }
      }
      w_gc = 0.5 * (lo + hi)
      let frac = (w_gc - prev_w) / (w - prev_w)
      let phase_at_gc = prev_phase + frac * (unwrapped_phase - prev_phase)
      pm_deg = 180.0 + phase_at_gc * 180.0 / @math.PI
    }

    // Check Phase Crossover: unwrapped_phase crosses -pi (-180 deg)
    let target_pi = -@math.PI
    if (prev_phase > target_pi && unwrapped_phase <= target_pi) ||
       (prev_phase < target_pi && unwrapped_phase >= target_pi) {
      let mut lo = prev_w
      let mut hi = w
      let mut lo_im = prev_val.im
      for iter = 0; iter < 50; iter = iter + 1 {
        let mid = 0.5 * (lo + hi)
        let mid_val = open_loop.eval_complex(@num.Complex::new(0.0, mid))
        if (lo_im >= 0.0 && mid_val.im >= 0.0) || (lo_im <= 0.0 && mid_val.im <= 0.0) {
          lo = mid
          lo_im = mid_val.im
        } else {
          hi = mid
        }
      }
      w_pc = 0.5 * (lo + hi)
      let val_pc = open_loop.eval_complex(@num.Complex::new(0.0, w_pc))
      let mag_pc = val_pc.abs()
      if mag_pc > 1.0e-15 {
        gm_db = -20.0 * @math.ln(mag_pc) / @math.ln(10.0)
      }
    }

    prev_w = w
    prev_val = val
    prev_mag = mag
    prev_phase = unwrapped_phase
  }

  {
    gain_margin_db: gm_db,
    phase_margin_deg: pm_deg,
    gain_crossover_freq: w_gc,
    phase_crossover_freq: w_pc,
    sensitivity_peak: max_ms,
  }
}

///|
pub fn StabilityMargins::to_string(self : StabilityMargins) -> String {
  "StabilityMargins(GM: \{self.gain_margin_db} dB @ \{self.phase_crossover_freq} rad/s, PM: \{self.phase_margin_deg} deg @ \{self.gain_crossover_freq} rad/s, Ms: \{self.sensitivity_peak})"
}

///|
pub impl Show for StabilityMargins with output(self, logger) {
  logger.write_string(self.to_string())
}
