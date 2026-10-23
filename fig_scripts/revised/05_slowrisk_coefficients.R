# ==============================================================================
# Main Figure 4: symptom-specific associations with continuous Slow Risk Load
# ==============================================================================

if (!exists("theme_revised")) {
  source("fig_scripts/revised/00_figure_setup.R")
}

ensure_analysis_data()

analysis_bin <- analysis_df %>%
  mutate(
    across(
      all_of(SYMPTOMS),
      ~ as.integer(.x > 0),
      .names = "{.col}_bin"
    )
  )

fit_one_symptom <- function(symptom_name, adjust_for_other_symptoms = FALSE) {
  response <- paste0(symptom_name, "_bin")

  predictors <- c(
    "slow_risk_z",
    "age",
    "gender_grp",
    "dutch_grp"
  )

  if (isTRUE(adjust_for_other_symptoms)) {
    predictors <- c(
      predictors,
      paste0(setdiff(SYMPTOMS, symptom_name), "_bin")
    )
  }

  formula <- reformulate(
    termlabels = predictors,
    response = response
  )

  model <- glm(
    formula,
    data = analysis_bin,
    family = binomial()
  )

  broom::tidy(model) %>%
    filter(term == "slow_risk_z") %>%
    transmute(
      symptom = symptom_name,
      model = ifelse(
        adjust_for_other_symptoms,
        "Covariates + other symptoms",
        "Covariates only"
      ),
      estimate = estimate,
      std_error = std.error,
      lower = estimate - stats::qnorm(0.975) * std.error,
      upper = estimate + stats::qnorm(0.975) * std.error,
      p_value = p.value
    )
}

coef_table <- bind_rows(
  purrr::map_dfr(
    SYMPTOMS,
    fit_one_symptom,
    adjust_for_other_symptoms = FALSE
  ),
  purrr::map_dfr(
    SYMPTOMS,
    fit_one_symptom,
    adjust_for_other_symptoms = TRUE
  )
)

symptom_order <- coef_table %>%
  filter(model == "Covariates + other symptoms") %>%
  arrange(estimate) %>%
  pull(symptom)

coef_table <- coef_table %>%
  mutate(
    symptom_label = unname(SYMPTOM_LABELS[symptom]),
    model = factor(
      model,
      levels = c(
        "Covariates only",
        "Covariates + other symptoms"
      )
    ),
    symptom_label = factor(
      symptom_label,
      levels = unname(SYMPTOM_LABELS[symptom_order])
    )
  )

write.csv(
  coef_table,
  "results/figure_data/figure_04_slowrisk_coefficients.csv",
  row.names = FALSE
)

x_min <- min(coef_table$lower, na.rm = TRUE)
x_max <- max(coef_table$upper, na.rm = TRUE)

# Add a small margin so the outer confidence intervals are not clipped
x_padding <- 0.04 * (x_max - x_min)

x_limits <- c(
  min(0, x_min - x_padding),
  x_max + x_padding
)

x_breaks <- pretty(
  x_limits,
  n = 5
)

fig_04 <- ggplot(
  coef_table,
  aes(
    x = estimate,
    y = symptom_label,
    color = model,
    shape = model
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.55,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = lower,
      xmax = upper
    ),
    height = 0.18,
    linewidth = 0.65,
    position = position_dodge(width = 0.52)
  ) +
  geom_point(
    size = 3.0,
    position = position_dodge(width = 0.52)
  ) +
  scale_color_manual(
    values = c(
      "Covariates only" = GREY_COLOR,
      "Covariates + other symptoms" = BLUE_COLOR
    )
  ) +
  scale_shape_manual(
    values = c(
      "Covariates only" = 17,
      "Covariates + other symptoms" = 16
    )
  ) +
  scale_x_continuous(
    breaks = x_breaks,
    limits = x_limits,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Symptom-specific associations with Slow Risk Load",
    x = paste(
      "Association with symptom presence",
      "(log-odds per 1 SD higher Slow Risk Load)",
      sep = "\n"
    ),
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    shape = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_revised +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.45
    )
  )

save_revised_figure(
  plot = fig_04,
  filename = "figure_04_slowrisk_coefficients",
  width = 9.5,
  height = 6.1
)
