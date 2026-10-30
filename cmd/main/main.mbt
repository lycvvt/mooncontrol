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
  let g1 = match @tf.TransferFunction::first_order(2.0, 0.5) {
    Ok(tf) => tf
    Err(e) => {
      println("Error creating G1: \{e}")
      return
    }
  }
  println("Plant G1(s):\{g1}")
  println("Poles of G1: \{format_poles(g1.poles())}")
  match g1.dc_gain() {
    Some(gain) => println("DC Gain: \{gain}")
    None => println("DC Gain: None")
  }
  println("")

  println("[Stage 2] Dimensionless Scaling & Asymptotic Step-Doubling Error Estimation")
  println("--------------------------------------------------------------------------------")
  // Plant A: Standard 2nd-order underdamped DC motor servo G(s) = 10 / (s^2 + 2s + 10)
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

  println("Step-Size Convergence & Asymptotic Step-Doubling Error Indicator (Plant A):")
  println("  dt (s) | eta=|lam|*dt | True Max Error     | Richardson Est     | Order p")
  println("  -------+--------------+--------------------+--------------------+---------")
  for r in study.records {
    let order_str = if r.order_estimate > 0.0 { "\{r.order_estimate}" } else { "---" }
    println("  \{r.dt}   | \{r.dimensionless_step} | \{r.max_error} | \{r.richardson_error_est} | \{order_str}")
  }
  println("  -------------------------------------------------------------------------")
  println("  Theoretical RK4 Order: \{study.theoretical_order} | Empirical Mean Order: \{study.mean_order}")

  // Demonstrate Adaptive RK45 / Step-Doubling Solver on a fast high-gain plant
  let fast_plant = @tf.TransferFunction::from_coeffs([100.0], [0.02, 1.0]).unwrap()
  let adap_res = @sim.simulate_step_adaptive(fast_plant, 0.5, atol=1.0e-9, rtol=1.0e-8).unwrap()
  println("  Adaptive Solver on Fast Plant G(s)=100/(0.02s+1) (lambda=-50, K=100):")
  println("    Accepted steps: \{adap_res.accepted_steps}, Rejected steps: \{adap_res.rejected_steps}, Max local error est: \{adap_res.max_estimated_local_error}")
  println("")

  println("[Stage 3] Exact 2D Stability Domain & Sampled-Data Closed-Loop Counterexamples")
  println("--------------------------------------------------------------------------------")
  let bounds_a = @sim.check_system_boundaries(plant_a)
  println("Plant A (Open-Loop Stable Underdamped Servo) Exact Stability Boundary:")
  println("  - Max Pole Magnitude |lambda|: \{bounds_a.max_pole_magnitude} rad/s")
  println("  - Directional Stability Radius: \{bounds_a.min_directional_radius} (along pole ray theta = 108.43 deg)")
  println("  - Exact Max Stable dt (RK4):   \{bounds_a.max_stable_dt} s")

  let plant_b = @tf.TransferFunction::from_coeffs([1.0], [1.0, 0.0, -9.8]).unwrap()
  let bounds_b = @sim.check_system_boundaries(plant_b)
  println("Plant B (Unstable Inverted Pendulum G(s) = 1/(s^2 - 9.8)):")
  println("  - Poles:                       \{format_poles(plant_b.poles())}")
  println("  - Open-Loop Hurwitz Stable:    \{bounds_b.is_open_loop_stable} -> Max Open-Loop Stable dt = \{bounds_b.max_stable_dt} s")

  let plant_1st = @tf.TransferFunction::from_coeffs([1.0], [1.0, 1.0]).unwrap()
  let rho_bad_dt = @sim.closed_loop_pid_spectral_radius(plant_1st, 8.0, 4.0, 0.2, 0.5).unwrap()
  let rho_good_dt = @sim.closed_loop_pid_spectral_radius(plant_1st, 8.0, 4.0, 0.2, 0.01).unwrap()
  let rho_bad_kd = @sim.closed_loop_pid_spectral_radius(plant_1st, 8.0, 4.0, 1.5, 0.01).unwrap()
  println("Sampled-Data Closed-Loop Transition Matrix Spectral Radius rho(Phi_CL(dt)):")
  println("  - Plant 1/(s+1), PID(8, 4, 0.2) at dt=0.50s (< open-loop 2.785s): rho = \{rho_bad_dt} (> 1: UNSTABLE!)")
  println("  - Plant 1/(s+1), PID(8, 4, 0.2) at dt=0.01s:                      rho = \{rho_good_dt} (< 1: STABLE)")
  println("  - Plant 1/(s+1), PID(8, 4, 1.5) at dt=0.01s (c*b*Kd = 1.5 > 1):   rho = \{rho_bad_kd} (> 1: ALGEBRAIC LOOP UNSTABLE!)")
  println("")

  println("[Stage 4] Quantitative Frequency Margins & State-Space Pole Placement")
  println("--------------------------------------------------------------------------------")
  let loop_tf = @tf.TransferFunction::from_coeffs([2.0], [1.0, 3.0, 3.0, 1.0]).unwrap()
  let margins = @analysis.compute_stability_margins(loop_tf)
  println("Open-Loop L(s) = 2 / (s + 1)^3 Quantitative Frequency Margins:")
  println("  - Gain Margin (GM):        \{margins.gain_margin_db} dB (at w_pc = \{margins.phase_crossover_freq} rad/s)")
  println("  - Phase Margin (PM):       \{margins.phase_margin_deg} deg (at w_gc = \{margins.gain_crossover_freq} rad/s)")
  println("  - Sensitivity Peak (M_s):  \{margins.sensitivity_peak}")

  let (a_pend, b_pend, _c_pend, ss_pend) = @analysis.analyze_transfer_function_state_space(plant_b).unwrap()
  let desired_poly = @num.Polynomial::new([1.0, 8.0, 32.0]) // target poles at -4 +- 4i
  let k_pend = @analysis.ackermann_pole_placement(a_pend, b_pend, desired_poly).unwrap()
  println("Inverted Pendulum State-Space Synthesis (G(s) = 1 / (s^2 - 9.8)):")
  println("  - Kalman Ranks:            Controllability = \{ss_pend.controllability_rank}/2, Observability = \{ss_pend.observability_rank}/2 (Minimal: \{ss_pend.is_minimal})")
  println("  - Ackermann Pole Placement to s = -4 +- 4i: State Feedback Gain K = [\{k_pend[0]}, \{k_pend[1]}]")
  println("")

  println("[Stage 5] Nonlinear & Delayed Closed-Loop Simulation (Beyond Symbolic Laplace)")
  println("--------------------------------------------------------------------------------")
  let pid_ctrl = @pid.PIDController::new(2.5, 4.0, 0.4, integral_limit=1.5, output_min=-1.5, output_max=1.5)
  println("Configured Discrete PID Controller: \{pid_ctrl}")
  let sim_no_aw = @sim.simulate_closed_loop_nonlinear(
    plant_a, 2.5, 4.0, 0.4, 5.0,
    setpoint=1.0, dt=0.005, u_min=-1.5, u_max=1.5,
    anti_windup_kb=0.0, transport_delay=0.02,
    disturbance_time=2.5, disturbance_magnitude=-0.4,
  ).unwrap()
  let sim_with_aw = @sim.simulate_closed_loop_nonlinear(
    plant_a, 2.5, 4.0, 0.4, 5.0,
    setpoint=1.0, dt=0.005, u_min=-1.5, u_max=1.5,
    anti_windup_kb=5.0, transport_delay=0.02,
    disturbance_time=2.5, disturbance_magnitude=-0.4,
  ).unwrap()
  let m_no_aw = @sim.measure_trajectory_metrics(sim_no_aw, setpoint=1.0)
  let m_with_aw = @sim.measure_trajectory_metrics(sim_with_aw, setpoint=1.0)
  println("Plant A with Saturation [-1.5, 1.5], Transport Delay tau=0.02s, Disturbance d=-0.4 at t=2.5s:")
  println("  - Without Back-Calc Anti-Windup (kb=0): Overshoot = \{m_no_aw.overshoot_percent}%, IAE = \{m_no_aw.iae}, ITAE = \{m_no_aw.itae}")
  println("  - With Back-Calc Anti-Windup    (kb=5): Overshoot = \{m_with_aw.overshoot_percent}%, IAE = \{m_with_aw.iae}, ITAE = \{m_with_aw.itae}")
  println("Ascii Nonlinear Anti-Windup Closed-Loop Response:")
  println(@plot.plot_ascii(sim_with_aw.times, sim_with_aw.outputs, width=50, height=10))

  println("[Stage 6] Zero-Dependency SVG Visualization Export")
  println("--------------------------------------------------------------------------------")
  let svg = @plot.plot_svg(
    "Nonlinear Anti-Windup PID Step & Disturbance Response",
    sim_with_aw.times,
    sim_with_aw.outputs,
    setpoint=1.0,
  )
  println("Generated SVG: \{svg.length()} bytes of pure scalable vector graphics.")
  println("All stages completed successfully!")
}
