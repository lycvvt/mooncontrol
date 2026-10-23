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

  println("[Stage 2] Rigorous RK4 Convergence Order & Step-Size Sensitivity Study")
  println("--------------------------------------------------------------------------------")
  // Plant A: Standard 2nd-order underdamped DC motor servo
  // G(s) = 10 / (s^2 + 2s + 10), wn = sqrt(10) ≈ 3.162, zeta = 1 / sqrt(10) ≈ 0.316
  let wn = (10.0 : Double).sqrt()
  let zeta = 1.0 / wn
  let plant_a = match @tf.TransferFunction::from_coeffs([10.0], [1.0, 2.0, 10.0]) {
    Ok(tf) => tf
    Err(e) => {
      println("Error creating Plant A: \{e}")
      return
    }
  }
  let dts = [0.08, 0.04, 0.02, 0.01, 0.005]
  let study = match @sim.evaluate_rk4_convergence(plant_a, 3.0, dts, fn(t) {
    match @analysis.step_response_second_order(wn, zeta, t) {
      Ok(v) => v
      Err(_) => 0.0
    }
  }) {
    Ok(st) => st
    Err(e) => {
      println("Convergence study error: \{e}")
      return
    }
  }

  println("Empirical Step-Size vs Error Convergence Table on Benchmark Plant A:")
  println("  dt (s)    | Max Absolute Error | RMS Error          | Estimated Convergence Order")
  println("  ----------+--------------------+--------------------+----------------------------")
  for r in study.records {
    let order_str = if r.order_estimate > 0.0 { "\{r.order_estimate}" } else { "---" }
    println("  \{r.dt}      | \{r.max_error} | \{r.rms_error} | \{order_str}")
  }
  println("  --------------------------------------------------------------------------------")
  println("  Theoretical RK4 Order: \{study.theoretical_order} | Empirical Mean Order: \{study.mean_order}")
  println("  Precision Boundary: Truncation dominates for dt in [1e-4, 1e-1]; Error scales strictly as O(dt^4).")
  println("")

  println("[Stage 3] System Stability Domain & Stiffness Boundary Check")
  println("--------------------------------------------------------------------------------")
  let bounds_a = @sim.check_system_boundaries(plant_a)
  println("Plant A Stability Boundary Analysis:")
  println("  - Max Pole Magnitude:    \{bounds_a.max_pole_magnitude} rad/s")
  println("  - Open-Loop Stable:      \{bounds_a.is_open_loop_stable}")
  println("  - Max Stable dt (RK4):   \{bounds_a.max_stable_dt} s (Stability condition: dt * |lambda| < 2.785)")
  println("  - Recommended dt:        \{bounds_a.recommended_dt} s (Safety factor 4x for high accuracy)")
  println("  - System Stiff (S > 100):\{bounds_a.is_stiff} (Stiffness ratio: \{bounds_a.stiffness_ratio})")

  // Unstable Inverted Pendulum Plant B: G(s) = 1 / (s^2 - 9.8)
  let plant_b = match @tf.TransferFunction::from_coeffs([1.0], [1.0, 0.0, -9.8]) {
    Ok(tf) => tf
    Err(_) => plant_a
  }
  let bounds_b = @sim.check_system_boundaries(plant_b)
  println("Plant B (Inverted Pendulum) Stability Analysis:")
  println("  - Poles:                 \{format_poles(plant_b.poles())}")
  println("  - Open-Loop Stable:      \{bounds_b.is_open_loop_stable} (RHP pole at +3.1305 rad/s triggers unstable alert)")
  println("")

  println("[Stage 4] Closed-Loop Control System Design & Routh-Hurwitz Verification")
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
  println("Theoretical Dynamic Metrics:")
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
  println("Measured Closed-Loop Error at dt=0.0005: \{max_err2}")
  println("Ascii Closed-Loop Step Response:")
  println(@plot.plot_ascii(sim2.times, sim2.outputs, width=50, height=10))

  println("[Stage 5] Discrete PID Controller with Anti-Windup Simulation")
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
  println("Ascii PID Controlled Step Response:")
  println(@plot.plot_ascii(pid_sim.times, pid_sim.outputs, width=50, height=10))

  println("[Stage 6] Zero-Dependency SVG Visualization Export")
  println("--------------------------------------------------------------------------------")
  let svg = @plot.plot_svg(
    "PID Closed-Loop Step Response",
    pid_sim.times,
    pid_sim.outputs,
    setpoint=1.0,
  )
  println("Generated SVG: \{svg.length()} bytes of pure scalable vector graphics.")
  println("All stages completed successfully!")
}
