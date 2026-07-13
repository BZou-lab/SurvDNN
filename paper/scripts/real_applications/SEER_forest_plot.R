# SEER forest plot

library(tidyverse)

# Import data
hcc_validdata <- readRDS("data/seer_validdata_unprocessed.rds")

# Stratify subjects into risk groups
risk_predictions <- readRDS("results/KM_SEER_71_selected.rds")
strata <- ifelse(
  risk_predictions[[1]] >= median(risk_predictions[[1]], na.rm = TRUE),
  "High",
  "Low"
)

# Add row for risk group to the test set
hcc_validdata$risk_group <- strata[rownames(hcc_validdata)]
if (all(is.na(hcc_validdata$risk_group)) && length(strata) == nrow(hcc_validdata)) {
  hcc_validdata$risk_group <- strata
}
hcc_validdata$risk_group <- factor(hcc_validdata$risk_group, levels = c("Low", "High"))

# Biomarker variables for forest plots ------------------------------------
# AFP_status (binary variable that is either "Positive" or "Negative")
# M_Stage (binary variable that is either M0 or M1)
# Radiation (binary variable that is either Yes or "No/Unknown")
# tumor_size (continuous variable for tumor size; units are mm)
# Age (continuous age in years)

forest_binary_vars <- c(
  "AFP_status",
  "M_Stage",
  "Radiation"
)

forest_continuous_vars <- c(
  "tumor_size",
  "Age"
)

forest_binary_event_levels <- c(
  AFP_status = "Positive",
  M_Stage = "M1",
  Radiation = "Yes"
)

forest_variable_labels <- c(
  AFP_status = "AFP Positive",
  Radiation = "Radiation Therapy",
  M_Stage = "M Stage",
  tumor_size = "Tumor Size",
  Age = "Age at Diagnosis"
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
  se <- sqrt((p_high * (1 - p_high) / n_high) + (p_low * (1 - p_low) / n_low))
  
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
  se <- sqrt((n_high + n_low) / (n_high * n_low) + smd^2 / (2 * (n_high + n_low - 2)))
  
  tibble(
    Variable = variable,
    SMD = smd,
    Lower = smd - 1.96 * se,
    Upper = smd + 1.96 * se,
    `Low-risk mean` = mean(low_values),
    `High-risk mean` = mean(high_values)
  )
}

seer_continuous_forest_data <- map_dfr(
  forest_continuous_vars,
  ~ compute_standardized_mean_difference(hcc_validdata, .x)
) %>%
  mutate(
    Variable_label = unname(forest_variable_labels[Variable]),
    Variable_label = fct_reorder(Variable_label, abs(SMD), .desc = FALSE)
  )

seer_binary_forest_data <- map_dfr(
  forest_binary_vars,
  ~ compute_prevalence_difference(hcc_validdata, .x)
) %>%
  mutate(
    Variable_label = unname(forest_variable_labels[Variable]),
    Variable_label = fct_reorder(Variable_label, abs(Difference), .desc = FALSE)
  )

binary_row_levels <- levels(seer_binary_forest_data$Variable_label)
radiation_index <- match("Radiation Therapy", binary_row_levels)
m_stage_index <- match("M Stage", binary_row_levels)
binary_row_levels[c(radiation_index, m_stage_index)] <-
  binary_row_levels[c(m_stage_index, radiation_index)]
seer_binary_forest_data$Variable_label <- factor(
  seer_binary_forest_data$Variable_label,
  levels = binary_row_levels
)

binary_direction_label_y <- nlevels(seer_binary_forest_data$Variable_label) + 0.9

continuous_x_limit <- max(abs(c(
  seer_continuous_forest_data$Lower,
  seer_continuous_forest_data$Upper
)), na.rm = TRUE)

binary_x_limit <- max(abs(c(
  seer_binary_forest_data$Lower,
  seer_binary_forest_data$Upper
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

seer_continuous_forest_plot <- ggplot(
  seer_continuous_forest_data,
  aes(x = SMD, y = Variable_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.05, linewidth = 0.6, color = "#009E73") +
  geom_point(size = 2.5, color = "#009E73") +
  scale_x_continuous(
    limits = c(-continuous_x_limit, continuous_x_limit),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    x = "Standardized mean difference"
  ) +
  forest_theme

seer_binary_forest_plot <- ggplot(
  seer_binary_forest_data,
  aes(x = Difference, y = Variable_label)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.05, linewidth = 0.6, color = "#009E73") +
  geom_point(size = 2.5, color = "#009E73") +
  scale_x_continuous(
    limits = c(-binary_x_limit, binary_x_limit),
    breaks = c(-30, -15, 0, 15, 30),
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
    x = -0.65 * binary_x_limit,
    xend = -0.75 * binary_x_limit,
    y = binary_direction_label_y,
    yend = binary_direction_label_y,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed"),
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
    x = 0.65 * binary_x_limit,
    xend = 0.75 * binary_x_limit,
    y = binary_direction_label_y,
    yend = binary_direction_label_y,
    arrow = grid::arrow(length = grid::unit(0.07, "inches"), type = "closed"),
    linewidth = 0.35,
    color = "#D55E00"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Prevalence difference, percentage points"
  ) +
  forest_theme +
  theme(plot.margin = ggplot2::margin(16, 10, 0, 8))

seer_binary_forest_grob <- ggplotGrob(seer_binary_forest_plot)
seer_continuous_forest_grob <- ggplotGrob(seer_continuous_forest_plot)
aligned_panel_widths <- grid::unit.pmax(
  seer_binary_forest_grob$widths,
  seer_continuous_forest_grob$widths
)
seer_binary_forest_grob$widths <- aligned_panel_widths
seer_continuous_forest_grob$widths <- aligned_panel_widths

panel_spacer_grob <- grid::rectGrob(
  gp = grid::gpar(col = NA, fill = NA)
)

seer_combined_forest_plot <- gridExtra::arrangeGrob(
  seer_binary_forest_grob,
  panel_spacer_grob,
  seer_continuous_forest_grob,
  ncol = 1,
  heights = c(length(forest_binary_vars), 0.32, length(forest_continuous_vars))
)

save_forest_plot <- function(plot, filename_stem, width = 7.5, height = 3.15) {
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

save_forest_plot(seer_combined_forest_plot, "SEER_biomarker_forest_plot")
