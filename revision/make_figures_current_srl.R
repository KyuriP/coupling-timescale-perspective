# ==============================================================================
# revision/make_figures_current_srl.R
#
# Re-creates the manuscript figures with the current Slow Risk Load and
# analytic sample (n = 23,270), using the existing figure scripts unchanged
# apart from their input/output paths (and axis labels where the plotted
# quantity changed).
#
#   Fig. 3   descriptives          archive/figure_scripts/legacy/fig_descriptives.R
#   Fig. 4   SRL effects on h      fig_scripts/revised/04_activation_differences.R
#   Supp S1  distributions         archive/figure_pipeline_legacy.R (section 14)
#   Supp S2  NCT                   fig_scripts/revised/03_nct_comparison.R
#   Supp S3  multigroup activation fig_scripts/revised/A1_activation_sensitivity.R
#
# Before running:
#   python3 revision/revision_analyses.py --sections C,mgtables,export,rexport
#   source("revision/nct_srl_split.R")
#
# Run from the repository root:
#
#   source("revision/make_figures_current_srl.R")
#
# ==============================================================================

source("R/00_setup.R")                          # theme_paper, group_colors, packages
source("fig_scripts/revised/00_figure_setup.R")  # theme_revised, save_revised_figure

dir.create("figs/revised", recursive = TRUE, showWarnings = FALSE)

run_script <- function(path, replacements = list(), from = NULL, to = NULL) {
  code <- readLines(path, warn = FALSE)
  if (!is.null(from)) {
    start <- grep(from, code, fixed = TRUE)[1]
    end <- if (is.null(to)) length(code) else grep(to, code, fixed = TRUE)[1]
    stopifnot(!is.na(start), !is.na(end))
    code <- code[start:end]
  }
  code <- paste(code, collapse = "\n")
  for (r in replacements) {
    n_hits <- lengths(regmatches(code, gregexpr(r[1], code, perl = TRUE)))
    if (n_hits == 0L) stop("Pattern not found in ", path, ": ", r[1])
    code <- gsub(r[1], r[2], code, perl = TRUE)
  }
  eval(parse(text = code), envir = globalenv())
}

# ------------------------------------------------------------------------------
# Analytic sample (written by revision_analyses.py --sections rexport)
# ------------------------------------------------------------------------------

analysis_df <- read.csv(file.path("data", "analysis_sample_current.csv")) %>%
  mutate(
    slow_group = factor(slow_group, levels = c("Low Slow Risk", "High Slow Risk")),
    dutch_grp = factor(dutch_grp, levels = c("Dutch", "Non-Dutch")),
    age_grp = factor(age_grp, levels = c("Younger", "Older"))
  )

# ------------------------------------------------------------------------------
# Fig. 3: descriptives
# ------------------------------------------------------------------------------

run_script("archive/figure_scripts/legacy/fig_descriptives.R")
ggsave("figs/revised/fig_descriptives_current.pdf", fig_descriptives,
       width = 7, height = 9, units = "in")
ggsave("figs/revised/fig_descriptives_current.png", fig_descriptives,
       width = 7, height = 9, units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# Supp. Fig. S1: distributions
# ------------------------------------------------------------------------------

run_script(
  "archive/figure_pipeline_legacy.R",
  from = "slow_dist_long <- analysis_df %>%",
  to = "plot_layout(heights = c(2, 1), axis_titles = \"collect\")"
)
ggsave("figs/revised/supp_distribution_hist_density_current.pdf", supp_dist2,
       width = 7, height = 7.5, units = "in")

# ------------------------------------------------------------------------------
# Fig. 4: SRL effects on symptom activation (moderated Ising model)
# ------------------------------------------------------------------------------

run_script(
  "fig_scripts/revised/04_activation_differences.R",
  list(
    c('"results",\\s*"network_invariance_exact",\\s*"slow_risk",\\s*"exact_multigroup_ising",\\s*"activation_differences_free_J_free_h_with_95CI\\.csv"',
      '"results", "revision_analyses", "figure_srl_activation_effects.csv"'),
    c("Higher in High Slow Risk", "Higher activation at higher SRL"),
    c("Higher in Low Slow Risk", "Lower activation at higher SRL"),
    c("Activation difference \\(High minus Low Slow Risk\\)",
      "Change in activation per SD of Slow Risk Load (log-odds)"),
    c("results/figure_data/figure_03_activation_differences\\.csv",
      "results/figure_data/figure_04_srl_activation_effects.csv"),
    c('filename = "figure_03_activation_differences"',
      'filename = "figure_04_srl_activation_effects"'),
    c("x_min <- min\\(activation\\$ci_95_lower\\) - 0\\.04",
      "x_min <- min(0, min(activation$ci_95_lower)) - 0.04")
  )
)

# ------------------------------------------------------------------------------
# Supp. Fig. S3: multigroup activation differences under two specifications
# ------------------------------------------------------------------------------

mg <- read.csv(file.path("results", "revision_analyses", "MG_activation_estimates.csv"))
sensitivity <- mg %>%
  transmute(
    symptom,
    model = dplyr::recode(model,
                          "fully group-specific" = "Group-specific interactions",
                          "common interactions" = "Common interactions"),
    difference = diff,
    lower = lo,
    upper = hi
  ) %>%
  mutate(
    symptom_label = unname(SYMPTOM_LABELS[symptom]),
    model = factor(model, levels = c("Group-specific interactions", "Common interactions"))
  )

run_script(
  "fig_scripts/revised/A1_activation_sensitivity.R",
  list(
    c("results/figure_data/supp_figure_S2_activation_sensitivity\\.csv",
      "results/figure_data/supp_figure_S3_activation_sensitivity_current.csv"),
    c("figs/figure_A1_activation_sensitivity",
      "figs/revised/supp_figure_S3_activation_sensitivity_current")
  ),
  from = "order_symptoms <- sensitivity %>%"
)

# ------------------------------------------------------------------------------
# Supp. Fig. S2: Network Comparison Test
# ------------------------------------------------------------------------------

nct_dir <- file.path("results", "revision_analyses", "nct")
if (file.exists(file.path(nct_dir, "nct_summary.csv"))) {
  # Summary aliases retained for compatibility with earlier figure templates.
  nct_sum <- read.csv(file.path(nct_dir, "nct_summary.csv"))
  nct_sum$M <- nct_sum$network_structure_statistic_M
  nct_sum$S <- nct_sum$global_strength_difference_S
  write.csv(nct_sum, file.path(nct_dir, "nct_summary_for_figure.csv"), row.names = FALSE)

  run_script(
    "fig_scripts/revised/03_nct_comparison.R",
    list(
      c("results/network_invariance_exact/slow_risk/nct_10000/nct_network_low_risk_10000\\.csv",
        "results/revision_analyses/nct/nct_network_low.csv"),
      c("results/network_invariance_exact/slow_risk/nct_10000/nct_network_high_risk_10000\\.csv",
        "results/revision_analyses/nct/nct_network_high.csv"),
      c("results/network_invariance_exact/slow_risk/nct_10000/nct_summary_10000\\.csv",
        "results/revision_analyses/nct/nct_summary_for_figure.csv"),
      c('"figs/revised/figure_nct_comparison\\.pdf"', '"figs/revised/figure_nct_comparison_current.pdf"'),
      c('"figs/revised/figure_nct_comparison\\.png"', '"figs/revised/figure_nct_comparison_current.png"')
    )
  )
  stopifnot(
    isTRUE(all.equal(M_value, nct_sum$M)),
    isTRUE(all.equal(S_value, nct_sum$S)),
    isTRUE(all.equal(GS_low, nct_sum$global_strength_group_1)),
    isTRUE(all.equal(GS_high, nct_sum$global_strength_group_2))
  )
} else {
  stop("NCT results are missing; run revision/nct_srl_split.R before making all figures.")
}

message(
  "\nFigures written to figs/revised/:\n",
  "  fig_descriptives_current.pdf\n",
  "  supp_distribution_hist_density_current.pdf\n",
  "  figure_04_srl_activation_effects.pdf\n",
  "  supp_figure_S3_activation_sensitivity_current.pdf\n",
  "  figure_nct_comparison_current.pdf\n"
)
