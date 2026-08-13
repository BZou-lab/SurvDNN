################ SEER Heatmap: Plot from final dataset ####################

library(ggplot2)
library(grid)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

make_seer_heatmap <- function(dfdraw) {
  
  ## Rename methods
  colnames(dfdraw)[colnames(dfdraw) == "XGB"] <- "XGBoost"
  colnames(dfdraw)[colnames(dfdraw) == "COX"] <- "Cox"
  
  ## Reorder methods
  method_order <- c("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")
  dfdraw <- dfdraw[, method_order, drop = FALSE]
  
  ## Rename features
  rownames(dfdraw)[rownames(dfdraw) == "Age"] <- "Age at Diagnosis"
  rownames(dfdraw)[rownames(dfdraw) == "Radiation"] <- "Radiation Therapy"
  rownames(dfdraw)[rownames(dfdraw) == "Single LN region"] <- "Single LN Region"
  rownames(dfdraw)[rownames(dfdraw) == "Surgery type"] <- "Surgery Type"
  rownames(dfdraw)[rownames(dfdraw) == "Time from diag to trt"] <- "Time to Treatment"
  rownames(dfdraw)[rownames(dfdraw) == "Tumor size"] <- "Tumor Size"
  rownames(dfdraw)[rownames(dfdraw) == "Grade"] <- "Tumor Grade"
  rownames(dfdraw)[rownames(dfdraw) == "Marital status"] <- "Marital Status"
  rownames(dfdraw)[rownames(dfdraw) == "First malignant primary"] <- "First Primary Tumor"
  rownames(dfdraw)[rownames(dfdraw) == "Number of malignant tumors"] <- "Number of Tumors"
  
  heat_df <- as.data.frame(as.table(dfdraw))
  colnames(heat_df) <- c("Feature", "Method", "Selected")
  
  heat_df$Feature <- factor(
    heat_df$Feature,
    levels = rev(rownames(dfdraw))
  )
  
  heat_df$Method <- factor(
    heat_df$Method,
    levels = method_order
  )
  
  ggplot(heat_df, aes(x = Method, y = Feature, fill = Selected)) +
    
    geom_tile(color = "white", linewidth = 0.35) +
    
    scale_fill_gradient(
      low = "#F7FBFF",
      high = "#08519C",
      limits = c(0, 100),
      name = "Times \nSelected",
    ) +
    
    labs(x = NULL, y = NULL) +
    
    theme_minimal(base_size = 12) +
    
    guides(
      fill = guide_colorbar(
        barwidth = unit(3, "mm"),
        barheight = unit(22, "mm")
      )
    ) +
    
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        size = 9,
        face = "bold",
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      
      axis.text.y = element_text(size = 9),
      
      legend.title = element_text(
        face = "bold",
        size = 10,
        margin = ggplot2::margin(b = 10)
      ),
      
      panel.border = element_rect(
        color = "grey80",
        fill = NA,
        linewidth = 0.5
      )
    )
}

dfdraw <- readRDS("results/SEER_Heatmap.rds")

p <- make_seer_heatmap(dfdraw)

p

ggsave(
  filename = "results/figures/SEER_Heatmap.pdf",
  plot = p,
  width = 4.25,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)


################ METABRIC Heatmap: Plot from final dataset ################

make_metabric_heatmap <- function(dfdraw) {
  
  ## Rename methods
  colnames(dfdraw)[colnames(dfdraw) == "XGB"] <- "XGBoost"
  colnames(dfdraw)[colnames(dfdraw) == "COX"] <- "Cox"
  
  ## Reorder methods
  method_order <- c("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")
  dfdraw <- dfdraw[, method_order, drop = FALSE]
  
  ## Rename features
  rownames(dfdraw)[rownames(dfdraw) == "Age"] <- 
    "Age at Diagnosis"
  
  rownames(dfdraw)[rownames(dfdraw) == "Positive Lymph Node #"] <- 
    "Positive Lymph Nodes"
  
  rownames(dfdraw)[rownames(dfdraw) == "Relapse"] <- 
    "Relapse Status"
  
  rownames(dfdraw)[rownames(dfdraw) == "BRCA Subtype"] <- 
    "BRCA Mutation"
  
  heat_df <- as.data.frame(as.table(dfdraw))
  colnames(heat_df) <- c("Feature", "Method", "Selected")
  
  heat_df$Feature <- factor(
    heat_df$Feature,
    levels = rev(rownames(dfdraw))
  )
  
  heat_df$Method <- factor(
    heat_df$Method,
    levels = method_order
  )
  
  ggplot(heat_df, aes(x = Method, y = Feature, fill = Selected)) +
    
    geom_tile(color = "white", linewidth = 0.35) +
    
    scale_fill_gradient(
      low = "#F7FBFF",
      high = "#08519C",
      limits = c(0, 100),
      name = "Times \nSelected",
    ) +
    
    labs(x = NULL, y = NULL) +
    
    theme_minimal(base_size = 12) +
    
    guides(
      fill = guide_colorbar(
        barwidth = unit(3, "mm"),
        barheight = unit(22, "mm")
      )
    ) +
    
    theme(
      panel.grid = element_blank(),
      
      axis.text.x = element_text(
        size = 9,
        face = "bold",
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      
      axis.text.y = element_text(size = 9),
      
      legend.title = element_text(
        face = "bold",
        size = 10,
        margin = ggplot2::margin(b = 10)
      ),
      
      panel.border = element_rect(
        color = "grey80",
        fill = NA,
        linewidth = 0.5
      )
    )
}

dfdraw <- readRDS("results/METABRIC_Heatmap.rds")

p <- make_metabric_heatmap(dfdraw)

p

ggsave(
  filename = "results/figures/METABRIC_Heatmap.pdf",
  plot = p,
  width = 4.25,
  height = 4.8,
  units = "in",
  device = cairo_pdf
)
