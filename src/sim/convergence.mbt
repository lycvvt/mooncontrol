///|
/// Convergence evaluation record for a single step size dt.
pub(all) struct StepConvergenceRecord {
  dt : Double
  max_error : Double
  rms_error : Double
  order_estimate : Double
} derive(Eq)

pub extend StepConvergenceRecord with Eq::{equal, not_equal}
pub extend StepConvergenceRecord with Show::{output}

///|
/// Comprehensive study of RK4 numerical convergence across varying step sizes.
pub(all) struct ConvergenceStudy {
  records : Array[StepConvergenceRecord]
  mean_order : Double
  theoretical_order : Double
  precision_floor : Double
} derive(Eq)

pub extend ConvergenceStudy with Eq::{equal, not_equal}
pub extend ConvergenceStudy with Show::{output}

///|
/// Run systematic step-size variation study to verify RK4 4th-order convergence O(dt^4).
/// dts: array of descending step sizes, e.g. [0.04, 0.02, 0.01, 0.005, 0.0025, 0.00125]
pub fn evaluate_rk4_convergence(
  tf : @tf.TransferFunction,
  duration : Double,
  dts : Array[Double],
  exact_fn : (Double) -> Double
) -> Result[ConvergenceStudy, String] {
  if dts.length() < 2 {
    return Err("At least two step sizes are required to evaluate convergence order")
  }

  let records = Array(capacity=dts.length())
  let mut prev_dt = 0.0
  let mut prev_err = 0.0
  let mut order_sum = 0.0
  let mut order_count = 0
  let mut min_err = 1.0e12

  for i = 0; i < dts.length(); i = i + 1 {
    let dt = dts[i]
    let sim = match simulate_step(tf, duration, dt~) {
      Ok(s) => s
      Err(e) => return Err(e)
    }

    let mut max_err = 0.0
    let mut sum_sq = 0.0
    let n = sim.times.length()

    for k = 0; k < n; k = k + 1 {
      let t = sim.times[k]
      let y_sim = sim.outputs[k]
      let y_exact = exact_fn(t)
      let diff = (y_sim - y_exact).abs()
      if diff > max_err {
        max_err = diff
      }
      sum_sq = sum_sq + diff * diff
    }

    let rms = (sum_sq / n.to_double()).sqrt()
    if max_err < min_err {
      min_err = max_err
    }

    let mut order = 0.0
    if i > 0 && max_err > 1.0e-13 && prev_err > 1.0e-13 {
      let log_err_ratio = @math.ln(prev_err / max_err)
      let log_dt_ratio = @math.ln(prev_dt / dt)
      order = log_err_ratio / log_dt_ratio
      order_sum = order_sum + order
      order_count = order_count + 1
    }

    records.push({ dt, max_error: max_err, rms_error: rms, order_estimate: order })
    prev_dt = dt
    prev_err = max_err
  }

  let mean_order = if order_count > 0 {
    order_sum / order_count.to_double()
  } else {
    4.0
  }

  Ok({
    records,
    mean_order,
    theoretical_order: 4.0,
    precision_floor: min_err,
  })
}

///|
pub fn StepConvergenceRecord::to_string(self : StepConvergenceRecord) -> String {
  "dt: \{self.dt}, max_err: \{self.max_error}, rms: \{self.rms_error}, order: \{self.order_estimate}"
}

///|
pub impl Show for StepConvergenceRecord with output(self, logger) {
  logger.write_string(self.to_string())
}

///|
pub fn ConvergenceStudy::to_string(self : ConvergenceStudy) -> String {
  let mut s = "ConvergenceStudy(theoretical_order: \{self.theoretical_order}, mean_order: \{self.mean_order}, floor: \{self.precision_floor}\n"
  for r in self.records {
    s = s + "  " + r.to_string() + "\n"
  }
  s + ")"
}

///|
pub impl Show for ConvergenceStudy with output(self, logger) {
  logger.write_string(self.to_string())
}
