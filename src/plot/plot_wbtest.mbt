test "plot ascii generates output with axes" {
  let times = [0.0, 1.0, 2.0, 3.0, 4.0]
  let values = [0.0, 0.63, 0.86, 0.95, 0.98]
  let ascii = plot_ascii(times, values, width=20, height=8)
  assert_true(ascii.length() > 0)
  assert_true(ascii.contains("*"))
  assert_true(ascii.contains("|"))
  assert_true(ascii.contains("+"))
}

test "plot svg generates valid svg xml tags" {
  let times = [0.0, 1.0, 2.0]
  let values = [0.0, 0.5, 1.0]
  let svg = plot_svg("Unit Step Response", times, values, setpoint=1.0)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("</svg>"))
  assert_true(svg.contains("Unit Step Response"))
  assert_true(svg.contains("<polyline"))
  assert_true(svg.contains("Setpoint (1)"))
}
