# ==============================================================================
# figure_nct_comparison.R
#
# NCT figure (Supplementary Fig. S2):
#   A. Low Slow Risk network
#   B. High Slow Risk network
#   C. Global strength summary
#
# Run from the repository root:
#   source("fig_scripts/revised/03_nct_comparison.R")
# The manuscript version is made by revision/make_figures_current_srl.R.
# ==============================================================================

# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------

required_packages <- c("qgraph")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Please install the following package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(qgraph)
})

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit) || length(hit) == 0L) {
    stop(
      "None of the candidate files were found:\n",
      paste("  -", paths, collapse = "\n")
    )
  }
  hit
}

read_network_matrix <- function(path) {
  x <- read.csv(path, row.names = 1, check.names = FALSE)
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  m
}

extract_first_numeric <- function(df, candidates) {
  for (nm in candidates) {
    if (nm %in% names(df)) {
      val <- suppressWarnings(as.numeric(df[[nm]][1]))
      if (!is.na(val)) return(val)
    }
  }
  NA_real_
}

# ------------------------------------------------------------------------------
# File locations
# ------------------------------------------------------------------------------

low_file <- find_first_existing(c(
  "results/network_invariance_exact/slow_risk/nct_10000/nct_network_low_risk_10000.csv",
  "results/network_invariance_exact/slow_risk/nct_network_low_risk_10000.csv",
  "nct_network_low_risk_10000.csv"
))

high_file <- find_first_existing(c(
  "results/network_invariance_exact/slow_risk/nct_10000/nct_network_high_risk_10000.csv",
  "results/network_invariance_exact/slow_risk/nct_network_high_risk_10000.csv",
  "nct_network_high_risk_10000.csv"
))

summary_file <- find_first_existing(c(
  "results/network_invariance_exact/slow_risk/nct_10000/nct_summary_10000.csv",
  "results/network_invariance_exact/slow_risk/nct_summary_10000.csv",
  "nct_summary_10000.csv"
))

dir.create("figs/revised", recursive = TRUE, showWarnings = FALSE)

pdf_out <- "figs/revised/figure_nct_comparison.pdf"
png_out <- "figs/revised/figure_nct_comparison.png"

# ------------------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------------------

network_low <- read_network_matrix(low_file)
network_high <- read_network_matrix(high_file)
summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)

if (!identical(dim(network_low), dim(network_high))) {
  stop("Low and High network matrices do not have the same dimensions.")
}

if (!identical(rownames(network_low), rownames(network_high))) {
  stop("Low and High network matrices do not have the same node ordering.")
}

# ------------------------------------------------------------------------------
# Labels and style
# ------------------------------------------------------------------------------

LOW_COLOR  <- "#1B9E77"
HIGH_COLOR <- "#D95F02"
DARK_GREY  <- "grey25"

abbr_map <- c(
  anh = "Anh",
  dep = "Dep",
  slp = "Slp",
  ene = "Ene",
  app = "App",
  glt = "Glt",
  con = "Con",
  mot = "Mot",
  sui = "Sui"
)

node_names <- rownames(network_low)
if (is.null(node_names)) {
  node_names <- c("anh","dep","slp","ene","app","glt","con","mot","sui")
}

node_labels <- abbr_map[node_names]
node_labels[is.na(node_labels)] <- node_names[is.na(node_labels)]

layout_shared <- qgraph::averageLayout(network_low, network_high)
maximum_edge <- max(abs(c(network_low, network_high)), na.rm = TRUE)

# ------------------------------------------------------------------------------
# Extract NCT summary values
# ------------------------------------------------------------------------------

M_value <- extract_first_numeric(summary_df, c(
  "network_structure_statistic_M",
  "network_structure_invariance",
  "structure_invariance",
  "M",
  "nwinv.real"
))

M_p <- extract_first_numeric(summary_df, c(
  "p_structure",
  "structure_invariance_p",
  "network_structure_p",
  "pval.structure",
  "pval.nwinv"
))

S_value <- extract_first_numeric(summary_df, c(
  "global_strength_difference_S",
  "global_strength_difference",
  "S",
  "glstrinv.real"
))

S_p <- extract_first_numeric(summary_df, c(
  "p_global_strength",
  "global_strength_p",
  "pval.glstrinv"
))

GS_low <- extract_first_numeric(summary_df, c(
  "global_strength_group_1",
  "global_strength_low",
  "glstr.low"
))

GS_high <- extract_first_numeric(summary_df, c(
  "global_strength_group_2",
  "global_strength_high",
  "glstr.high"
))

# stop if anything is missing from the summary file
if (any(!is.finite(c(GS_low, GS_high, M_value, M_p, S_value, S_p)))) {
  stop("NCT summary is missing a required finite statistic or p-value.")
}

# ------------------------------------------------------------------------------
# Plotting function
# ------------------------------------------------------------------------------

draw_nct_figure <- function() {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
         heights = c(1.0, 0.78))
  
  # --------------------------------------------------------------------------
  # Panel A: Low Slow Risk network
  # --------------------------------------------------------------------------
  par(mar = c(1.2, 1.2, 3.2, 1.2), family = "")
  qgraph::qgraph(
    network_low,
    layout = layout_shared,
    labels = node_labels,
    label.cex = 1.15,
    edge.color = LOW_COLOR,
    border.color = LOW_COLOR,
    border.width = 2,
    vsize = 7.8,
    esize = 8,
    maximum = maximum_edge,
    minimum = 0,
    cut = 0,
    details = FALSE,
    legend = FALSE,
    title = paste0("A  Low Slow Risk\nGlobal strength = ", sprintf("%.2f", GS_low))
  )
  
  # --------------------------------------------------------------------------
  # Panel B: High Slow Risk network
  # --------------------------------------------------------------------------
  par(mar = c(1.2, 1.2, 3.2, 1.2), family = "")
  qgraph::qgraph(
    network_high,
    layout = layout_shared,
    labels = node_labels,
    label.cex = 1.15,
    edge.color = HIGH_COLOR,
    border.color = HIGH_COLOR,
    border.width = 2,
    vsize = 7.8,
    esize = 8,
    maximum = maximum_edge,
    minimum = 0,
    cut = 0,
    details = FALSE,
    legend = FALSE,
    title = paste0("B  High Slow Risk\nGlobal strength = ", sprintf("%.2f", GS_high))
  )
  
  # --------------------------------------------------------------------------
  # Panel C: Global strength summary
  # --------------------------------------------------------------------------
  par(mar = c(4.5, 3.2, 3.4, 2.0), family = "")
  
  x_min <- 21
  x_max <- 25
  x_mid <- mean(c(GS_low, GS_high))
  
  # Put dots essentially on the x-axis
  y_axis_bottom <- 0.52
  y_dumbbell <- 0.48
  
  # Put labels and values above the dots
  y_label <- 0.95
  y_number <- 0.84
  
  plot(
    NA,
    xlim = c(x_min, x_max),
    ylim = c(y_axis_bottom, 1.45),
    xlab = "Global strength",
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    bty = "n",
    main = "C  Global strength"
  )
  
  # Draw the axis first
  axis(
    side = 1,
    at = seq(x_min, x_max, by = 1),
    labels = seq(x_min, x_max, by = 1),
    tck = -0.025
  )
  
  axis_y <- par("usr")[3]
  y_dumbbell <- axis_y
  
  # Dashed vertical guides
  segments(
    x0 = c(GS_low, GS_high),
    y0 = y_number - 0.03,
    x1 = c(GS_low, GS_high),
    y1 = y_dumbbell,
    col = "grey75",
    lty = 3,
    lwd = 1,
    xpd = NA
  )
  
  # Horizontal dumbbell connector
  segments(
    x0 = GS_low,
    y0 = y_dumbbell,
    x1 = GS_high,
    y1 = y_dumbbell,
    lwd = 2.2,
    col = "grey55",
    xpd = NA
  )
  
  # Dots drawn last so they appear in front of the axis and connector
  points(
    x = c(GS_low, GS_high),
    y = rep(y_dumbbell, 2),
    pch = 21,
    bg = c(LOW_COLOR, HIGH_COLOR),
    col = "black",
    cex = 2,
    lwd = 1,
    xpd = NA
  )
  
  # Group labels
  text(
    GS_low,
    y_label,
    labels = "Low Slow Risk",
    col = LOW_COLOR,
    cex = 0.92
  )
  
  text(
    GS_high,
    y_label,
    labels = "High Slow Risk",
    col = HIGH_COLOR,
    cex = 0.92
  )
  
  # Exact values directly below the labels
  text(
    c(GS_low, GS_high),
    y_number,
    labels = sprintf("%.2f", c(GS_low, GS_high)),
    col = c(LOW_COLOR, HIGH_COLOR),
    font = 2,
    cex = 0.92
  )
  
  # NCT statistics
  text(
    x_mid,
    1.34,
    labels = sprintf(
      "Structure invariance: M = %.3f, p = %.3f",
      M_value,
      M_p
    ),
    cex = 0.95,
    col = DARK_GREY
  )
  
  text(
    x_mid,
    1.23,
    labels = sprintf(
      "Global strength difference: S = %.3f, p = %.4f",
      S_value,
      S_p
    ),
    cex = 0.95,
    col = DARK_GREY
  )
}

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

grDevices::cairo_pdf(pdf_out, width = 9.4, height = 7.6)
draw_nct_figure()
grDevices::dev.off()

grDevices::png(
  filename = png_out,
  width = 3000,
  height = 2400,
  res = 320,
  bg = "white"
)
draw_nct_figure()
grDevices::dev.off()

message("Saved: ", pdf_out)
message("Saved: ", png_out)
