
library(ggplot2)
library(dplyr)
library(gridExtra)

summary_s2 = readRDS("results/SimulationS1_result_CIndex.rds")
summary_s3 = readRDS("results/SimulationS2_result_CIndex.rds")

# ============================================================================
# Colors
# ============================================================================

base_method_levels <- c("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")

method_levels <- c(
  "SurvDeepFIT",
  "Cox",
  "AFT",
  "RSF",
  "XGBoost"
)

method_colors <- c(
  "SurvDeepFIT" = "#D55E00",  
  "Cox"     = "#56B4E9",  
  "AFT"     = "#009E73",  
  "RSF"     = "#E69F00",  
  "XGBoost" = "#CC79A7"   
)

method_labels <- c(
  "SurvDeepFIT" = "SurvDeepFIT",
  "Cox" = "Cox",
  "AFT" = "AFT",
  "RSF" = "RSF",
  "XGBoost" = "XGBoost"
)

# ============================================================================
# Custom labeller using expression for rho
# sim_label: "S1" or "S2"
# ============================================================================

# (facet labels use Unicode \u03C1 for rho)

# ============================================================================
# Plot function
# sim_label: "S1" or "S2" (used in title and facet strips)
# ============================================================================

make_row_plot <- function(data, sample_size, sim_label, show_legend = FALSE) {
  
  plot_data <- data %>%
    filter(n == sample_size, Type == "Selected") %>%
    mutate(
      Method = factor(Method, levels = base_method_levels),
      Plot_Method = paste(as.character(Method)),
      Plot_Method = factor(Plot_Method, levels = method_levels)
    )
  
  # Create facet label with Unicode rho
  plot_data$facet_label <- paste0("Scenario ", sim_label, ": n = ", sample_size,
                                  ", \u03C1 = ", plot_data$rho)
  plot_data$facet_label <- factor(
    plot_data$facet_label,
    levels = paste0("Scenario ", sim_label, ": n = ", sample_size,
                    ", \u03C1 = ", c(0, 0.3, 0.5))
  )
  
  p <- ggplot(plot_data, aes(x = Method, y = mean_CIndex, fill = Plot_Method)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(
      aes(ymin = mean_CIndex - sd_CIndex, ymax = mean_CIndex + sd_CIndex),
      width = 0.25, linewidth = 0.4
    ) +
    facet_wrap(~ facet_label, nrow = 1) +
    scale_fill_manual(values = method_colors, breaks = method_levels, labels = method_labels) +
    coord_cartesian(ylim = c(0.5, 0.9)) +
    labs(x = NULL, y = "C-index", fill = "Method") +
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_text(size = 15),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "grey85"),
      strip.text = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12, face = "plain",
                                margin = margin(t = 2, r = 0, b = 2, l = 0)),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    ) +
    ggtitle(bquote("Scenario " * .(sim_label) * " (n = " * .(sample_size) * ")"))
  
  if (show_legend) {
    p <- p + theme(legend.position = "bottom", legend.direction = "horizontal")
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

get_legend <- function(plot) {
  plot_grob <- ggplotGrob(plot)
  legend_index <- which(sapply(plot_grob$grobs, function(grob) grob$name) == "guide-box")
  plot_grob$grobs[[legend_index]]
}

# ============================================================================
# Create and save plots
# ============================================================================

# --- s2 plot (Simulation S1) ---
p_s2_1000 <- make_row_plot(summary_s2, 1000, "1", show_legend = FALSE)
p_s2_3000 <- make_row_plot(summary_s2, 3000, "1", show_legend = FALSE)

# --- s3 plot (Simulation S2) ---
p_s3_1000 <- make_row_plot(summary_s3, 1000, "2", show_legend = FALSE)
p_s3_3000 <- make_row_plot(summary_s3, 3000, "2", show_legend = FALSE) +
  theme(plot.margin = margin(t = 5.5, r = 5.5, b = 18, l = 5.5))
plot_legend <- get_legend(make_row_plot(summary_s3, 3000, "2", show_legend = TRUE))

# --- combined plot (Scenario 1 on top, Scenario 2 below; one legend at bottom) ---
plot_cindex <- grid.arrange(
  p_s2_1000,
  p_s2_3000,
  p_s3_1000,
  p_s3_3000,
  plot_legend,
  nrow = 5,
  heights = c(1, 1, 1, 1, 0.25)
)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = "results/figures/Simulation_CIndex.pdf",
  plot = plot_cindex,
  width = 9,
  height = 10.5,
  units = "in",
  device = cairo_pdf
)
