///|
/// Example 2: Second-Order Plant with Closed-Loop Discrete PID Control
fn main {
  println("=======================================================")
  println("   Example 2: Second-Order System & PID Control        ")
  println("=======================================================")

  // Open-loop plant: G(s) = 2.0 / (s^2 + 1.5s + 0.5)
  let plant = match @tf.TransferFunction::from_coeffs([2.0], [1.0, 1.5, 0.5]) {
    Ok(tf) => tf
    Err(e) => {
      println("Failed to create plant: \{e}")
      return
    }
  }

  println("Open-Loop Plant G(s):\{plant}")

  // Open-loop step response simulation
  let open_sim = match @sim.simulate_step(plant, 10.0, dt=0.001) {
    Ok(s) => s
    Err(e) => {
      println("Simulation failed: \{e}")
      return
    }
  }
  println("Open-Loop Steady-State Value: \{open_sim.outputs[open_sim.outputs.length() - 1]}")
  println("\nOpen-Loop Step Curve:")
  println(@plot.plot_ascii(open_sim.times, open_sim.outputs, width=45, height=8))

  // Design PID Controller: Kp = 3.0, Ki = 1.2, Kd = 0.8 with anti-windup clamping to 2.5
  let pid = @pid.PIDController::new(3.0, 1.2, 0.8, integral_limit=2.5)
  println("Configured Controller: \{pid}")

  // Run closed loop simulation with target setpoint = 1.0
  let closed_sim = match @sim.simulate_closed_loop_pid(plant, pid, 8.0, setpoint=1.0, dt=0.001) {
    Ok(s) => s
    Err(e) => {
      println("Closed-loop simulation failed: \{e}")
      return
    }
  }

  println("\nClosed-Loop PID Step Response (Setpoint = 1.0):")
  println(@plot.plot_ascii(closed_sim.times, closed_sim.outputs, width=45, height=12))

  // Export SVG
  let svg = @plot.plot_svg("Closed-Loop PID Tracking Response", closed_sim.times, closed_sim.outputs, setpoint=1.0)
  println("Generated SVG: \{svg.length()} bytes.")
}
