# ==============================================================================
# Main Figure 4: Slow Risk Load effects on activation and interactions
#
# Coefficients of continuous SRL (per SD) in the exact moderated Ising model in
# which SRL moderates all activation and interaction parameters.
#   a  activation of the nine symptoms
#   b  the 36 pairwise interactions
# Both panels are on the same log-odds scale. Orange marks effects that remain
# after Holm correction. Bars are 95% Wald intervals.
#
# Input: results/revision_analyses/C_srl_h.csv
#        results/revision_analyses/C_srl_J.csv
#        (python3 revision/revision_analyses.py --sections C)
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_figure_setup.R")
}

h_file <- file.path("results", "revision_analyses", "C_srl_h.csv")
J_file <- file.path("results", "revision_analyses", "C_srl_J.csv")

require_result_file(h_file)
require_result_file(J_file)

activation <- read.csv(h_file, stringsAsFactors = FALSE) %>%
  transmute(
    label = unname(SYMPTOM_LABELS[symptom]),
    est, lo, hi,
    holm = p_holm < 0.05
  ) %>%
  arrange(est) %>%
  mutate(label = factor(label, levels = label))

interactions <- read.csv(J_file, stringsAsFactors = FALSE) %>%
  mutate(
    v1 = sub("--.*", "", edge),
    v2 = sub(".*--", "", edge)
  ) %>%
  transmute(
    label = paste0(SYMPTOM_LABELS[v1], " – ", tolower(SYMPTOM_LABELS[v2])),
    est, lo, hi,
    holm = p_holm < 0.05
  ) %>%
  arrange(est) %>%
  mutate(label = factor(label, levels = label))

write.csv(
  bind_rows(
    mutate(activation, parameter = "activation"),
    mutate(interactions, parameter = "interaction")
  ),
  "results/figure_data/figure_04_srl_activation_effects.csv",
  row.names = FALSE
)

x_limits <- c(
  min(0, activation$lo, interactions$lo) - 0.03,
  max(activation$hi, interactions$hi) + 0.03
)

effect_panel <- function(dat, title) {
  ggplot(
    dat,
    aes(
      x = est,
      y = label,
      color = holm,
      fill = holm
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.6,
      color = "grey45"
    ) +
    geom_segment(
      aes(
        x = lo,
        xend = hi,
        yend = label
      ),
      linewidth = 0.8,
      alpha = 0.55
    ) +
    geom_point(
      shape = 21,
      size = 2.6,
      stroke = 0.8
    ) +
    scale_color_manual(
      values = c(`TRUE` = HIGH_COLOR, `FALSE` = "grey55"),
      guide = "none"
    ) +
    scale_fill_manual(
      values = c(`TRUE` = scales::alpha(HIGH_COLOR, 0.35), `FALSE` = "white"),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq(-0.2, 0.8, by = 0.2)
    ) +
    coord_cartesian(
      xlim = x_limits,
      clip = "off"
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme_revised +
    theme(
      plot.title = element_text(face = "bold", size = 12.5, hjust = 0),
      plot.title.position = "plot",
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.45),
      axis.text.y = element_text(size = 9.5)
    )
}

panel_a <- effect_panel(
  activation,
  "a  Activation of the nine symptoms"
)

panel_b <- effect_panel(
  interactions,
  "b  Interactions between the 36 symptom pairs"
) +
  labs(
    x = "Change per SD of Slow Risk Load (log-odds)"
  )

# Heights proportional to the number of rows, so both panels share the same
# row spacing
fig_04 <- panel_a / panel_b +
  plot_layout(heights = c(nrow(activation), nrow(interactions)))

save_revised_figure(
  plot = fig_04,
  filename = "figure_04_srl_activation_effects",
  width = 7.5,
  height = 10.5
)
