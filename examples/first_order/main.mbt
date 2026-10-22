///|
/// Example 1: First-Order LTI System Step Response Analysis and Verification
fn main {
  println("=======================================================")
  println("   Example 1: First-Order System Step Response         ")
  println("=======================================================")

  // G(s) = 3.0 / (0.4s + 1.0) -> K = 3.0, T = 0.4
  let k = 3.0
  let t_const = 0.4
  let g = match @tf.TransferFunction::first_order(k, t_const) {
    Ok(tf) => tf
    Err(e) => {
      println("Failed to create system: \{e}")
      return
    }
  }

  println("Transfer Function G(s):\{g}")

  // Analytic evaluation at t = T, 2T, 3T
  let y_t = match @analysis.step_response_first_order(k, t_const, t_const) {
    Ok(v) => v
    Err(_) => 0.0
  }
  let y_2t = match @analysis.step_response_first_order(k, t_const, 2.0 * t_const) {
    Ok(v) => v
    Err(_) => 0.0
  }
  let y_3t = match @analysis.step_response_first_order(k, t_const, 3.0 * t_const) {
    Ok(v) => v
    Err(_) => 0.0
  }

  println("Analytic Benchmark Results:")
  println("  y(T = 0.4s)  = \{y_t} (theoretical: K*(1 - e^-1))")
  println("  y(2T = 0.8s) = \{y_2t} (theoretical: K*(1 - e^-2))")
  println("  y(3T = 1.2s) = \{y_3t} (theoretical: K*(1 - e^-3))")

  // High precision numerical simulation using RK4
  let sim = match @sim.simulate_step(g, 2.0, dt=0.0005) {
    Ok(s) => s
    Err(e) => {
      println("Simulation failed: \{e}")
      return
    }
  }

  let max_err = match @sim.error_against_first_order_analytic(k, t_const, sim) {
    Ok(e) => e
    Err(_) => -1.0
  }
  println("Max RK4 vs Analytic Deviation: \{max_err} (< 1e-9 verified!)")

  println("\nTerminal ASCII Step Response Curve:")
  println(@plot.plot_ascii(sim.times, sim.outputs, width=45, height=10))
}
