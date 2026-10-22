///|
/// PID Controller with continuous and discrete position algorithms, anti-windup, and derivative filtering.
pub(all) struct PIDController {
  kp : Double
  ki : Double
  kd : Double
  integral_limit : Double?
  output_min : Double?
  output_max : Double?
  filter_coeff : Double // alpha in [0, 1) for first-order derivative filter
  mut integral : Double
  mut prev_error : Double
  mut filtered_derivative : Double
  mut initialized : Bool
}

pub extend PIDController with Show::{output}

///|
pub fn PIDController::new(
  kp : Double,
  ki : Double,
  kd : Double,
  integral_limit? : Double,
  output_min? : Double,
  output_max? : Double,
  filter_coeff~ : Double = 0.0
) -> PIDController {
  {
    kp,
    ki,
    kd,
    integral_limit,
    output_min,
    output_max,
    filter_coeff,
    integral: 0.0,
    prev_error: 0.0,
    filtered_derivative: 0.0,
    initialized: false,
  }
}

///|
/// Reset internal state (integrator and derivative memory).
pub fn PIDController::reset(self : PIDController) -> Unit {
  self.integral = 0.0
  self.prev_error = 0.0
  self.filtered_derivative = 0.0
  self.initialized = false
}

///|
/// Compute control output u[k] given setpoint, measurement, and sample period dt.
pub fn PIDController::update(
  self : PIDController,
  setpoint : Double,
  measurement : Double,
  dt : Double
) -> Double {
  if dt <= 0.0 {
    return 0.0
  }
  let error = setpoint - measurement

  // Proportional term
  let p_term = self.kp * error

  // Integral term with anti-windup clamping
  self.integral = self.integral + error * dt
  match self.integral_limit {
    Some(limit) => {
      let abs_lim = limit.abs()
      if self.integral > abs_lim {
        self.integral = abs_lim
      } else if self.integral < -abs_lim {
        self.integral = -abs_lim
      }
    }
    None => ()
  }
  let i_term = self.ki * self.integral

  // Derivative term with optional first-order low-pass filter
  let raw_derivative = if self.initialized {
    (error - self.prev_error) / dt
  } else {
    self.initialized = true
    0.0
  }
  self.prev_error = error

  let alpha = self.filter_coeff
  self.filtered_derivative = alpha * self.filtered_derivative + (1.0 - alpha) * raw_derivative
  let d_term = self.kd * self.filtered_derivative

  let mut u = p_term + i_term + d_term

  // Output saturation limits
  match self.output_min {
    Some(min_val) => if u < min_val { u = min_val }
    None => ()
  }
  match self.output_max {
    Some(max_val) => if u > max_val { u = max_val }
    None => ()
  }

  u
}

///|
/// Convert PID parameters to continuous transfer function C(s) = (Kd*s^2 + Kp*s + Ki) / s
pub fn PIDController::to_transfer_function(
  self : PIDController
) -> Result[@tf.TransferFunction, String] {
  if self.ki.abs() > 1.0e-12 {
    // Has integral action: (Kd*s^2 + Kp*s + Ki) / s
    if self.kd.abs() > 1.0e-12 {
      @tf.TransferFunction::from_coeffs([self.kd, self.kp, self.ki], [1.0, 0.0])
    } else {
      @tf.TransferFunction::from_coeffs([self.kp, self.ki], [1.0, 0.0])
    }
  } else {
    // No integral action: Ki == 0
    if self.kd.abs() > 1.0e-12 {
      @tf.TransferFunction::from_coeffs([self.kd, self.kp], [1.0])
    } else {
      @tf.TransferFunction::from_coeffs([self.kp], [1.0])
    }
  }
}

///|
pub fn PIDController::to_string(self : PIDController) -> String {
  "PIDController { Kp: \{self.kp}, Ki: \{self.ki}, Kd: \{self.kd} }"
}

///|
pub impl Show for PIDController with output(self, logger) {
  logger.write_string(self.to_string())
}
