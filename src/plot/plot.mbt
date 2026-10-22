///|
/// Render a lightweight ASCII terminal graph from time series data.
pub fn plot_ascii(
  times : Array[Double],
  values : Array[Double],
  width~ : Int = 60,
  height~ : Int = 15
) -> String {
  let len = times.length()
  if len == 0 || values.length() == 0 {
    return "No data to plot"
  }

  let mut min_y = values[0]
  let mut max_y = values[0]
  for v in values {
    if v < min_y { min_y = v }
    if v > max_y { max_y = v }
  }

  // Ensure non-zero range
  if (max_y - min_y).abs() < 1.0e-9 {
    max_y = min_y + 1.0
  }

  let t_max = times[len - 1]
  let grid = Array::make(height, Array::make(width, ' '))
  for row = 0; row < height; row = row + 1 {
    grid[row] = Array::make(width, ' ')
  }

  for i = 0; i < len; i = i + 1 {
    let t = times[i]
    let y = values[i]
    let col = ((t / t_max) * (width - 1).to_double()).to_int()
    let norm_y = (y - min_y) / (max_y - min_y)
    let row = ((1.0 - norm_y) * (height - 1).to_double()).to_int()

    let clamped_col = if col < 0 { 0 } else if col >= width { width - 1 } else { col }
    let clamped_row = if row < 0 { 0 } else if row >= height { height - 1 } else { row }
    grid[clamped_row][clamped_col] = '*'
  }

  let buf = StringBuilder()
  buf.write_string(" ^\n")
  for row = 0; row < height; row = row + 1 {
    let y_val = max_y - (row.to_double() / (height - 1).to_double()) * (max_y - min_y)
    let y_str = "\{y_val}"
    // Truncate to 6 chars
    let padded = if y_str.length() > 6 {
      y_str[0:6].to_owned()
    } else {
      let pad = String::make(6 - y_str.length(), ' ')
      "\{pad}\{y_str}"
    }
    buf.write_string("\{padded} |")
    for col = 0; col < width; col = col + 1 {
      buf.write_char(grid[row][col])
    }
    buf.write_char('\n')
  }
  buf.write_string("       +")
  for col = 0; col < width; col = col + 1 {
    buf.write_char('-')
  }
  buf.write_string(">\n")
  buf.write_string("       0.0")
  let t_str = "\{t_max}s"
  let spaces = width - 3 - t_str.length()
  if spaces > 0 {
    buf.write_string(String::make(spaces, ' '))
  }
  buf.write_string("\{t_str}\n")

  buf.to_string()
}

///|
/// Generate crisp, standalone SVG document for step response curves.
pub fn plot_svg(
  title : String,
  times : Array[Double],
  values : Array[Double],
  setpoint? : Double,
  width~ : Int = 800,
  height~ : Int = 500
) -> String {
  let len = times.length()
  if len == 0 || values.length() == 0 {
    return "<svg xmlns=\"http://www.w3.org/2000/svg\"><text>No data</text></svg>"
  }

  let mut min_y = values[0]
  let mut max_y = values[0]
  for v in values {
    if v < min_y { min_y = v }
    if v > max_y { max_y = v }
  }
  match setpoint {
    Some(sp) => {
      if sp < min_y { min_y = sp }
      if sp > max_y { max_y = sp }
    }
    None => ()
  }

  // Margin padding
  let pad_margin = (max_y - min_y) * 0.1
  let y_lower = min_y - pad_margin
  let y_upper = max_y + pad_margin
  let y_span = if (y_upper - y_lower).abs() < 1.0e-9 { 1.0 } else { y_upper - y_lower }
  let t_max = if times[len - 1] > 0.0 { times[len - 1] } else { 1.0 }

  let margin_left = 70.0
  let margin_right = 40.0
  let margin_top = 60.0
  let margin_bottom = 50.0

  let plot_w = width.to_double() - margin_left - margin_right
  let plot_h = height.to_double() - margin_top - margin_bottom

  let buf = StringBuilder()
  buf.write_string(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 \{width} \{height}\" width=\"\{width}\" height=\"\{height}\" style=\"background-color:#ffffff; font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;\">\n",
  )

  // Title
  buf.write_string(
    "  <text x=\"\{width / 2}\" y=\"35\" font-size=\"18\" font-weight=\"bold\" text-anchor=\"middle\" fill=\"#1e293b\">\{title}</text>\n",
  )

  // Plot border and background
  buf.write_string(
    "  <rect x=\"\{margin_left}\" y=\"\{margin_top}\" width=\"\{plot_w}\" height=\"\{plot_h}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1\" />\n",
  )

  // Horizontal grid lines (5 subdivisions)
  for g = 0; g <= 5; g = g + 1 {
    let ratio = g.to_double() / 5.0
    let gy = margin_top + plot_h * (1.0 - ratio)
    let gval = y_lower + ratio * y_span
    buf.write_string(
      "  <line x1=\"\{margin_left}\" y1=\"\{gy}\" x2=\"\{margin_left + plot_w}\" y2=\"\{gy}\" stroke=\"#e2e8f0\" stroke-dasharray=\"4,4\" />\n",
    )
    buf.write_string(
      "  <text x=\"\{margin_left - 10.0}\" y=\"\{gy + 4.0}\" font-size=\"11\" text-anchor=\"end\" fill=\"#64748b\">\{gval}</text>\n",
    )
  }

  // Vertical grid lines (5 subdivisions)
  for g = 0; g <= 5; g = g + 1 {
    let ratio = g.to_double() / 5.0
    let gx = margin_left + plot_w * ratio
    let gval = ratio * t_max
    buf.write_string(
      "  <line x1=\"\{gx}\" y1=\"\{margin_top}\" x2=\"\{gx}\" y2=\"\{margin_top + plot_h}\" stroke=\"#e2e8f0\" stroke-dasharray=\"4,4\" />\n",
    )
    buf.write_string(
      "  <text x=\"\{gx}\" y=\"\{margin_top + plot_h + 20.0}\" font-size=\"11\" text-anchor=\"middle\" fill=\"#64748b\">\{gval}s</text>\n",
    )
  }

  // Setpoint line
  match setpoint {
    Some(sp) => {
      let sp_y = margin_top + plot_h * (1.0 - (sp - y_lower) / y_span)
      buf.write_string(
        "  <!-- Setpoint line -->\n  <line x1=\"\{margin_left}\" y1=\"\{sp_y}\" x2=\"\{margin_left + plot_w}\" y2=\"\{sp_y}\" stroke=\"#ef4444\" stroke-width=\"1.5\" stroke-dasharray=\"6,3\" />\n",
      )
      buf.write_string(
        "  <text x=\"\{margin_left + plot_w - 5.0}\" y=\"\{sp_y - 6.0}\" font-size=\"11\" fill=\"#ef4444\" text-anchor=\"end\">Setpoint (\{sp})</text>\n",
      )
    }
    None => ()
  }

  // Response polyline points
  buf.write_string("  <!-- Response Curve -->\n  <polyline fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.5\" stroke-linejoin=\"round\" stroke-linecap=\"round\" points=\"")
  for i = 0; i < len; i = i + 1 {
    let t = times[i]
    let y = values[i]
    let px = margin_left + (t / t_max) * plot_w
    let py = margin_top + (1.0 - (y - y_lower) / y_span) * plot_h
    if i > 0 {
      buf.write_char(' ')
    }
    buf.write_string("\{px},\{py}")
  }
  buf.write_string("\" />\n")

  // Axis Labels
  buf.write_string(
    "  <text x=\"\{margin_left + plot_w / 2.0}\" y=\"\{height.to_double() - 10.0}\" font-size=\"13\" font-weight=\"500\" text-anchor=\"middle\" fill=\"#334155\">Time (seconds)</text>\n",
  )
  buf.write_string(
    "  <text transform=\"rotate(-90 20 \{margin_top + plot_h / 2.0})\" x=\"20\" y=\"\{margin_top + plot_h / 2.0}\" font-size=\"13\" font-weight=\"500\" text-anchor=\"middle\" fill=\"#334155\">Amplitude</text>\n",
  )

  buf.write_string("</svg>\n")
  buf.to_string()
}
