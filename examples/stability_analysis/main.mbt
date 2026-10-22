///|
/// Example 3: Routh-Hurwitz and Pole-Placement Stability Analysis
fn main {
  println("=======================================================")
  println("   Example 3: Comprehensive Stability Analysis         ")
  println("=======================================================")

  // System A (Stable): s^3 + 6s^2 + 11s + 6 = 0
  let p_stable = @num.Polynomial::new([1.0, 6.0, 11.0, 6.0])
  let tf_stable = match @tf.TransferFunction::from_coeffs([1.0], [1.0, 6.0, 11.0, 6.0]) {
    Ok(tf) => tf
    Err(_) => return
  }

  println("\n[System A] Characteristic Polynomial: \{p_stable}")
  let rh_a = match @analysis.routh_hurwitz(p_stable) {
    Ok(rh) => rh
    Err(e) => {
      println("Routh analysis failed: \{e}")
      return
    }
  }
  println("Routh Array:\n\{rh_a}")
  println("Poles:")
  for pole in tf_stable.poles() {
    println("  pole: \{pole}")
  }
  println("Is Stable (Routh-Hurwitz): \{rh_a.is_stable}")
  println("Is Stable (Pole Real Part): \{@analysis.is_stable_by_poles(tf_stable)}")

  // System B (Unstable): s^3 + s^2 + 2s + 24 = 0
  let p_unstable = @num.Polynomial::new([1.0, 1.0, 2.0, 24.0])
  let tf_unstable = match @tf.TransferFunction::from_coeffs([1.0], [1.0, 1.0, 2.0, 24.0]) {
    Ok(tf) => tf
    Err(_) => return
  }

  println("\n[System B] Characteristic Polynomial: \{p_unstable}")
  let rh_b = match @analysis.routh_hurwitz(p_unstable) {
    Ok(rh) => rh
    Err(e) => {
      println("Routh analysis failed: \{e}")
      return
    }
  }
  println("Routh Array:\n\{rh_b}")
  println("Poles:")
  for pole in tf_unstable.poles() {
    println("  pole: \{pole}")
  }
  println("Is Stable (Routh-Hurwitz): \{rh_b.is_stable}")
  println("Is Stable (Pole Real Part): \{@analysis.is_stable_by_poles(tf_unstable)}")
  println("Unstable RHP Poles Count: \{rh_b.sign_changes}")
}
