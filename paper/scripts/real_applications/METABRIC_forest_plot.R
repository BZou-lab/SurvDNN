# METABRIC forest plot

library(tidyverse)

# Import data
metabric_validdata <- readRDS("data/validdata_45837.rds")

# Stratify subjects into risk groups
risk_predictions <- readRDS("results/KM_METABRIC_19_selected.rds")
strata <- ifelse(
  risk_predictions[[1]] >= median(risk_predictions[[1]], na.rm = TRUE),
  "High",
  "Low"
)

# Add row for risk group to the test set
metabric_validdata$risk_group <- strata[rownames(metabric_validdata)]
if (all(is.na(metabric_validdata$risk_group)) && length(strata) == nrow(metabric_validdata)) {
  metabric_validdata$risk_group <- strata
}
metabric_validdata$risk_group <- factor(metabric_validdata$risk_group, levels = c("Low", "High"))

# Biomarker variables for forest plots ------------------------------------
# ER_STATUS_P (0/1 binary)
# TP53 (0/1 binary)
# PR_STATUS_P (0/1 binary)
# NPI (continuous)
# TUMOR_SIZE (continuous)

forest_binary_vars <- c(
  "ER_STATUS_P",
  "TP53",
  "PR_STATUS_P"
)

forest_continuous_vars <- c(
  "NPI",
  "TUMOR_SIZE"
)

forest_binary_event_levels <- c(
  ER_STATUS_P = "1",
  TP53 = "1",
  PR_STATUS_P = "1"
)

forest_variable_labels <- c(
  ER_STATUS_P = "ER Positive",
  TP53 = "TP53 Mutation",
  PR_STATUS_P = "PR Positive",
  NPI = "NPI",
  TUMOR_SIZE = "Tumor Size"
)

compute_prevalence_difference <- function(data, variable) {
  event_level <- unname(forest_binary_event_levels[variable])
  
  low_values <- data %>%
    filter(risk_group == "Low") %>%
    pull(all_of(variable)) %>%
    as.character()
  
  high_values <- data %>%
    filter(risk_group == "High") %>%
    pull(all_of(variable)) %>%
    as.character()
  
  n_low <- sum(!is.na(low_values))
  n_high <- sum(!is.na(high_values))
  p_low <- mean(low_values == event_level, na.rm = TRUE)
  p_high <- mean(high_values == event_level, na.rm = TRUE)
  difference <- 100 * (p_high - p_low)
  se <- sqrt(
    (p_high * (1 - p_high) / n_high) +
      (p_low * (1 - p_low) / n_low)
  )
  
  tibble(
    Variable = variable,
    Difference = difference,
    Lower = difference - 1.96 * 100 * se,
    Upper = difference + 1.96 * 100 * se,
    `Low-risk prevalence` = 100 * p_low,
    `High-risk prevalence` = 100 * p_high
  )
}

compute_standardized_mean_difference <- function(data, variable) {
  low_values <- data %>%
    filter(risk_group == "Low") %>%
    pull(all_of(variable)) %>%
    as.numeric()
  
  high_values <- data %>%
    filter(risk_group == "High") %>%
    pull(all_of(variable)) %>%
    as.numeric()
  
  low_values <- low_values[!is.na(low_values)]
  high_values <- high_values[!is.na(high_values)]
  n_low <- length(low_values)
  n_high <- length(high_values)
  pooled_sd <- sqrt(
    ((n_high - 1) * var(high_values) + (n_low - 1) * var(low_values)) /
      (n_high + n_low - 2)
  )
  smd <- (mean(high_values) - mean(low_values)) / pooled_sd
  se <- sqrt(
    (n_high + n_low) / (n_high * n_low) +
      smd^2 / (2 * (n_high + n_low - 2))
  )
  
  tibble(
    Variable = variable,
    SMD = smd,
    Lower = smd - 1.96 * se,
    Upper = smd + 1.96 * se,
    `Low-risk mean` = mean(low_values),
    `High-risk mean` = mean(high_values)
  )
}

metabric_binary_forest_data <- map_dfr(
  forest_binary_vars,
  ~ compute_prevalence_difference(metabric_validdata, .x)
) %>%
  mutate(
    Variable_label = unname(forest_variable_labels[Variable]),
    Variable_label = fct_reorder(Variable_label, abs(Difference), .desc = FALSE)
  )

metabric_continuous_forest_data <- map_dfr(
  forest_continuous_vars,
  ~ compute_standardized_mean_difference(metabric_validdata, .x)
) %>%
  mutate(
    Variable_label = unname(forest_variable_labels[Variable]),
    Variable_label = factor(
      Variable_label,
      levels = c("Tumor Size", "NPI")
    )
  )

binary_direction_label_y <-
  nlevels(metabric_binary_forest_data$Variable_label) + 0.9

binary_x_limit <- max(abs(c(
  metabric_binary_forest_data$Lower,
  metabric_binary_forest_data$Upper
)), na.rm = TRUE)

continuous_x_limit <- max(abs(c(
  metabric_continuous_forest_data$Lower,
  metabric_continuous_forest_data$Upper
)), na.rm = TRUE)

forest_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black", size = 9.5),
    axis.title.x = element_text(
      color = "black",
      size = 10,
      margin = ggplot2::margin(t = 5)
    ),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    plot.margin = ggplot2::margin(0, 10, 0, 8)
  )

metabric_binary_forest_plot <- ggplot(
  metabric_binary_forest_data,
  aes(x = Difference, y = Variable_label)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray45",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper),
    height = 0.05,
    linewidth = 0.6,
    color = "#009E73"
  ) +
  geom_point(size = 2.5, color = "#009E73") +
  scale_x_continuous(
    limits = c(-25, 25),
    breaks = c(-20, -10, 0, 10, 20),
    labels = c("-20", "-10", "0", "10", "20"),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  annotate(
    "text",
    x = 0,
    y = binary_direction_label_y,
    label = "Enriched in Low-Risk Group",
    hjust = 1.12,
    size = 3.1,
    color = "#0072B2"
  ) +
  annotate(
    "segment",
    x = -16,
    xend = -18.5,
    y = binary_direction_label_y,
    yend = binary_direction_label_y,
    arrow = grid::arrow(
      length = grid::unit(0.07, "inches"),
      type = "closed"
    ),
    linewidth = 0.35,
    color = "#0072B2"
  ) +
  annotate(
    "text",
    x = 0,
    y = binary_direction_label_y,
    label = "Enriched in High-Risk Group",
    hjust = -0.12,
    size = 3.1,
    color = "#D55E00"
  ) +
  annotate(
    "segment",
    x = 16,
    xend = 18.5,
    y = binary_direction_label_y,
    yend = binary_direction_label_y,
    arrow = grid::arrow(
      length = grid::unit(0.07, "inches"),
      type = "closed"
    ),
    linewidth = 0.35,
    color = "#D55E00"
  ) +
  coord_cartesian(clip = "off") +
  labs(x = "Prevalence difference, percentage points") +
  forest_theme +
  theme(plot.margin = ggplot2::margin(16, 10, 0, 8))

metabric_continuous_forest_plot <- ggplot(
  metabric_continuous_forest_data,
  aes(x = SMD, y = Variable_label)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray45",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = Lower, xmax = Upper),
    height = 0.05,
    linewidth = 0.6,
    color = "#009E73"
  ) +
  geom_point(size = 2.5, color = "#009E73") +
  scale_x_continuous(
    limits = c(-continuous_x_limit, continuous_x_limit),
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(x = "Standardized mean difference") +
  forest_theme

# Align the plotting panels so the zero lines share the same horizontal center.
metabric_binary_forest_grob <- ggplotGrob(metabric_binary_forest_plot)
metabric_continuous_forest_grob <- ggplotGrob(metabric_continuous_forest_plot)
aligned_panel_widths <- grid::unit.pmax(
  metabric_binary_forest_grob$widths,
  metabric_continuous_forest_grob$widths
)
metabric_binary_forest_grob$widths <- aligned_panel_widths
metabric_continuous_forest_grob$widths <- aligned_panel_widths

panel_spacer_grob <- grid::rectGrob(
  gp = grid::gpar(col = NA, fill = NA)
)

metabric_combined_forest_plot <- gridExtra::arrangeGrob(
  metabric_binary_forest_grob,
  panel_spacer_grob,
  metabric_continuous_forest_grob,
  ncol = 1,
  heights = c(
    length(forest_binary_vars),
    0.32,
    length(forest_continuous_vars)
  )
)

save_forest_plot <- function(
  plot,
  filename_stem,
  width = 7.5,
  height = 3.15
) {
  dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
  
  ggsave(
    filename = paste0("results/figures/", filename_stem, ".pdf"),
    plot = plot,
    width = width,
    height = height
  )
  
  ggsave(
    filename = paste0("results/figures/", filename_stem, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = 600
  )
}

save_forest_plot(
  metabric_combined_forest_plot,
  "METABRIC_biomarker_forest_plot"
)

