################ METABRIC C-index and IBS: Plot from final dataset ################

library(ggplot2)
library(dplyr)
library(grid)
library(svglite)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# C-index --------------------------------------------------------------

# ============================================================================
# Colors and ordering
# ============================================================================

base_method_levels <- c("SurvDNN", "Cox", "AFT", "RSF", "XGBoost")

method_levels <- c(
  "SurvDNN",
  "Cox",
  "AFT",
  "RSF",
  "XGBoost"
)

method_colors <- c(
  "SurvDNN" = "#D55E00",  
  "Cox"     = "#56B4E9",  
  "AFT"     = "#009E73",  
  "RSF"     = "#E69F00",  
  "XGBoost" = "#CC79A7"   
)

method_labels <- c(
  "SurvDNN"  = "SurvDNN",
  "Cox"      = "Cox",
  "AFT"      = "AFT",
  "RSF"      = "RSF",
  "XGBoost"  = "XGBoost"
)

format_prediction_results <- function(df) {
  df %>%
    mutate(
      Method = recode(base, "Xgb" = "XGBoost"),
      Method = factor(Method, levels = base_method_levels),
      Type = if_else(grepl("-Selected$", Name), "Selected", "Original"),
      Type = factor(Type, levels = c("Original", "Selected"))
    ) %>%
    filter(Type == "Selected") %>%
    mutate(
      Plot_Method = paste(as.character(Method)),
      Plot_Method = factor(Plot_Method, levels = method_levels)
    )
}

make_prediction_plot <- function(df, y_label, y_limits) {
  bar_w <- 0.6
  
  ggplot(df, aes(x = Method, y = C_Index, fill = Plot_Method)) +
    geom_col(
      width = bar_w
    ) +
    geom_errorbar(
      aes(ymin = C_Index - SD, ymax = C_Index + SD),
      width = 0.15
    ) +
    scale_fill_manual(
      values = method_colors,
      breaks = method_levels,
      labels = method_labels,
      name = "Method"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = FALSE)) +
    coord_cartesian(ylim = y_limits) +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = c(0.35, 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.key.size = unit(3, "mm"),
      legend.spacing.x = unit(0.75, "mm"),
      legend.spacing.y = unit(0.5, "mm"),
      axis.title.y = element_text(size = 15),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "plain")
    ) +
    ggtitle("")
}

# ============================================================================
# C-index
# ============================================================================

df <- readRDS("results/METABRIC_CIndex.rds")
df <- format_prediction_results(df)

P3 <- make_prediction_plot(
  df = df,
  y_label = "C-index",
  y_limits = c(0.5, 0.85)
)

P3

ggsave(
  filename = "results/figures/METABRIC_CIndex.pdf",
  plot = P3,
  width = 4.55,
  height = 4.6,
  units = "in",
  device = cairo_pdf
)

# IBS -------------------------------------------------------------------------

# ============================================================================
# Colors and ordering
# ============================================================================

base_method_levels <- c("SurvDNN", "Cox", "AFT", "RSF", "XGBoost")

method_levels <- c(
  "SurvDNN",
  "Cox",
  "AFT",
  "RSF",
  "XGBoost"
)

method_colors <- c(
  "SurvDNN" = "#D55E00",  
  "Cox"     = "#56B4E9",  
  "AFT"     = "#009E73",  
  "RSF"     = "#E69F00",  
  "XGBoost" = "#CC79A7"   
)

method_labels <- c(
  "SurvDNN"  = "SurvDNN",
  "Cox"      = "Cox",
  "AFT"      = "AFT",
  "RSF"      = "RSF",
  "XGBoost"  = "XGBoost"
)

format_prediction_results <- function(df) {
  df %>%
    mutate(
      Method = recode(base, "Xgb" = "XGBoost"),
      Method = factor(Method, levels = base_method_levels),
      Type = if_else(grepl("-Selected$", Name), "Selected", "Original"),
      Type = factor(Type, levels = c("Original", "Selected"))
    ) %>%
    filter(Type == "Selected") %>%
    mutate(
      Plot_Method = paste(as.character(Method)),
      Plot_Method = factor(Plot_Method, levels = method_levels)
    )
}

make_prediction_plot <- function(df, y_label, y_limits) {
  bar_w <- 0.6
  
  ggplot(df, aes(x = Method, y = C_Index, fill = Plot_Method)) +
    geom_col(
      width = bar_w
    ) +
    geom_errorbar(
      aes(ymin = C_Index - SD, ymax = C_Index + SD),
      width = 0.15
    ) +
    scale_fill_manual(
      values = method_colors,
      breaks = method_levels,
      labels = method_labels,
      name = "Method"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = FALSE)) +
    coord_cartesian(ylim = y_limits) +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = c(0.35, 0.5),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 8.5),
      legend.key.size = unit(3, "mm"),
      legend.spacing.x = unit(0.75, "mm"),
      legend.spacing.y = unit(0.5, "mm"),
      axis.title.y = element_text(size = 15),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "plain")
    ) +
    ggtitle("")
}

# ============================================================================
# Integrated Brier Score (IBS)
# ============================================================================

df <- readRDS("results/METABRIC_IBS.rds")
df <- format_prediction_results(df)

P3 <- make_prediction_plot(
  df = df,
  y_label = "IBS",
  y_limits = c(0.0, 0.2)
)

P3

ggsave(
  filename = "results/figures/METABRIC_IBS.pdf",
  plot = P3,
  width = 4.6,
  height = 4.6,
  units = "in",
  device = cairo_pdf
)
