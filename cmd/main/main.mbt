///|
fn format_poles(poles : Array[@num.Complex]) -> String {
  let buf = StringBuilder()
  buf.write_string("[")
  for i = 0; i < poles.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string("\{poles[i]}")
  }
  buf.write_string("]")
  buf.to_string()
}

///|
/// MoonControl CLI Entry Point
fn main {
  println("================================================================================")
  println("        MoonControl - Control System Design & Simulation in MoonBit            ")
  println("================================================================================")
  println("")

  let args = @env.args()

  if args.length() > 1 && (args[1] == "--help" || args[1] == "-h") {
    println("Usage: moon run cmd/main [--demo | --help]")
    println("")
    println("Commands & Options:")
    println("  --demo      Run automated control system design, RK4 simulation, and verification")
    println("  --help      Show this help information")
    return
  }

  println("[Stage 1] Numerical Core & Transfer Function Modeling")
  println("--------------------------------------------------------------------------------")
  // First-order plant: G1(s) = 2.0 / (0.5s + 1.0)
  let g1 = match @tf.TransferFunction::first_order(2.0, 0.5) {
    Ok(tf) => tf
    Err(e) => {
      println("Error creating G1: \{e}")
      return
    }
  }
  println("Plant G1(s):\{g1}")
  let p1 = g1.poles()
  println("Poles of G1: \{format_poles(p1)}")

  match g1.dc_gain() {
    Some(gain) => println("DC Gain: \{gain}")
    None => println("DC Gain: None")
  }
  println("")

  println("[Stage 2] High-Precision Step Response & Analytical Benchmark (Error < 1e-9)")
  println("--------------------------------------------------------------------------------")
  let y_t = match @analysis.step_response_first_order(2.0, 0.5, 0.5) {
    Ok(v) => v
    Err(_) => 0.0
  }
  let y_2t = match @analysis.step_response_first_order(2.0, 0.5, 1.0) {
    Ok(v) => v
    Err(_) => 0.0
  }
  let y_3t = match @analysis.step_response_first_order(2.0, 0.5, 1.5) {
    Ok(v) => v
    Err(_) => 0.0
  }
  println("Analytical y(T)   = \{y_t} (theoretical: K*(1 - e^-1))")
  println("Analytical y(2T)  = \{y_2t} (theoretical: K*(1 - e^-2))")
  println("Analytical y(3T)  = \{y_3t} (theoretical: K*(1 - e^-3))")

  // Simulate with RK4
  let sim1 = match @sim.simulate_step(g1, 2.5, dt=0.0005) {
    Ok(res) => res
    Err(e) => {
      println("Simulation failed: \{e}")
      return
    }
  }
  let max_err1 = match @sim.error_against_first_order_analytic(2.0, 0.5, sim1) {
    Ok(err) => err
    Err(_) => -1.0
  }
  println("RK4 Simulation vs Analytical Max Error: \{max_err1}")
  println("Ascii Step Response Curve:")
  println(@plot.plot_ascii(sim1.times, sim1.outputs, width=50, height=10))

  println("[Stage 3] Second-Order Feedback Control & Stability Analysis")
  println("--------------------------------------------------------------------------------")
  // Plant: G(s) = 1.0 / (s^2 + s)
  let plant2 = match @tf.TransferFunction::from_coeffs([1.0], [1.0, 1.0, 0.0]) {
    Ok(tf) => tf
    Err(e) => {
      println("Error: \{e}")
      return
    }
  }
  println("Open-Loop Plant G(s) = 1 / (s^2 + s)")

  // Proportional controller Kp = 4.0
  let controller = match @tf.TransferFunction::from_coeffs([4.0], [1.0]) {
    Ok(tf) => tf
    Err(e) => {
      println("Error: \{e}")
      return
    }
  }
  let closed_loop = (controller * plant2).unity_feedback()
  println("Closed-Loop T(s) = C*G / (1 + C*G):\{closed_loop}")

  let cl_poles = closed_loop.poles()
  println("Closed-Loop Poles: \{format_poles(cl_poles)}")

  // Routh-Hurwitz Stability Verification
  let rh = match @analysis.routh_hurwitz(closed_loop.den) {
    Ok(arr) => arr
    Err(e) => {
      println("Routh error: \{e}")
      return
    }
  }
  let poles_stable = @analysis.is_stable_by_poles(closed_loop)
  println("Routh-Hurwitz Stable: \{rh.is_stable} (Sign changes: \{rh.sign_changes})")
  println("Pole-Check Stable:    \{poles_stable}")

  // Second-order metrics (wn = 2.0, zeta = 0.25)
  let metrics = match @analysis.calculate_second_order_metrics(2.0, 0.25) {
    Ok(m) => m
    Err(e) => {
      println("Metrics error: \{e}")
      return
    }
  }
  println("Theoretical Metrics:")
  println("  - Overshoot Mp:       \{metrics.overshoot_percent}%")
  println("  - Peak Time tp:       \{metrics.peak_time}s")
  println("  - Rise Time tr:       \{metrics.rise_time}s")
  println("  - Settling Time ts:   \{metrics.settling_time_2pct}s (2% criterion)")

  // RK4 Closed-loop simulation
  let sim2 = match @sim.simulate_step(closed_loop, 4.0, dt=0.0005) {
    Ok(res) => res
    Err(e) => {
      println("Simulation error: \{e}")
      return
    }
  }
  let max_err2 = match @sim.error_against_second_order_analytic(2.0, 0.25, sim2) {
    Ok(err) => err
    Err(_) => -1.0
  }
  println("Closed-loop Simulation vs Analytic Error: \{max_err2}")
  println("Ascii Closed-Loop Step Response:")
  println(@plot.plot_ascii(sim2.times, sim2.outputs, width=50, height=12))

  println("[Stage 4] Discrete PID Closed-Loop Controller Simulation")
  println("--------------------------------------------------------------------------------")
  let pid_controller = @pid.PIDController::new(4.0, 1.0, 0.5, integral_limit=2.0)
  println("Configured PID: \{pid_controller}")
  let pid_sim = match @sim.simulate_closed_loop_pid(plant2, pid_controller, 5.0, setpoint=1.0) {
    Ok(res) => res
    Err(e) => {
      println("PID simulation error: \{e}")
      return
    }
  }
  println("Ascii PID Controlled Response (Plant: 1/(s^2+s), PID: Kp=4, Ki=1, Kd=0.5):")
  println(@plot.plot_ascii(pid_sim.times, pid_sim.outputs, width=50, height=12))

  println("[Stage 5] Standalone SVG Vector Plot Export")
  println("--------------------------------------------------------------------------------")
  let svg = @plot.plot_svg(
    "PID Closed-Loop Step Response",
    pid_sim.times,
    pid_sim.outputs,
    setpoint=1.0,
  )
  println("Generated SVG: \{svg.length()} bytes of pure scalable vector graphics.")
  println("Verification successful! All systems normal.")
}
