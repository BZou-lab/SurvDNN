
library(survival)
library(tidyverse)
library(simsurv)
library(MASS)
library(deepTL)
library(survivalmodels)
library(randomForestSRC)
library(survivalsvm)
library(reticulate)
library(Hmisc)
library(xgboost)
library(survivalsvm)
library(pheatmap)
library(timeROC)
library(survminer)
library(gridExtra)

############ SEER ###############

pred_org = readRDS("results/KM_SEER_71_org.rds")

pred_selected = readRDS("results/KM_SEER_71_selected.rds")


args1 = 71
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) *2399

p_thershold = 0.05

use_seed = 256 +start.seed
set.seed(use_seed)

trainN = paste("data/traindata","_",use_seed,".rds",sep = "")
validN = paste("data/validdata","_",use_seed,".rds",sep = "")

Ptrain = readRDS(trainN)
Pvalid = readRDS(validN)
SEER_data_train = importDnnetSurv(x = Ptrain[,c(1:20,23:27)],
                                  y = Ptrain[,"OS_MONTHS"],
                                  e = Ptrain[,"OS_STATUS"])
SEER_data_valid = importDnnetSurv(x = Pvalid[,c(1:20,23:27)],
                                  y = Pvalid[,"OS_MONTHS"],
                                  e = Pvalid[,"OS_STATUS"])

### Risk stratification (using SurvDeepFIT with selected features only)
data_selected = data.frame(y = SEER_data_valid@y,
                      e = SEER_data_valid@e,
                      strata = ifelse(pred_selected[[1]] >= median(pred_selected[[1]]),"High","Low")) %>%
  mutate(strata = factor(strata, levels = c("Low", "High")))
fit_km = survfit(Surv(y,e)~strata, data = data_selected)

p1 <- ggsurvplot(
  fit_km,
  data = data_selected,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = FALSE,
  censor = FALSE,
  size = 1,
  palette = c("#0072B2", "#D55E00"),
  legend.title = "Risk Group",
  legend.labs = c("Low", "High"),
  xlab = "Time (months)",
  ylab = "Survival probability",
  break.time.by = 10,
  ggtheme = theme_classic(base_size = 11)
)$plot +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1))) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_blank()
  )

### Histological grade (3-tier scale) -----------------------------------------------------------
data_km <- data.frame(
  y = SEER_data_valid@y,
  e = SEER_data_valid@e,
  grade_2 = SEER_data_valid@x[, "grade_2"],
  grade_3 = SEER_data_valid@x[, "grade_3"]
) %>%
  mutate(
    strata = case_when(
      grade_2 == 0 & grade_3 == 0 ~ "I",
      grade_2 == 1 & grade_3 == 0 ~ "II",
      grade_2 == 0 & grade_3 == 1 ~ "III"
    ),
    strata = factor(strata, levels = c("I", "II", "III"))
  ) %>%
  filter(!is.na(strata))

fit_km <- survfit(Surv(y, e) ~ strata, data = data_km)

p2 <- ggsurvplot(
  fit_km,
  data = data_km,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = FALSE,
  censor = FALSE,
  size = 1,
  palette = c("#0072B2", "#D55E00", "#009E73"),
  legend.title = "Tumor Grade",
  legend.labs = c("I", "II", "III"),
  xlab = "Time (months)",
  ylab = "Survival probability",
  break.time.by = 10,
  ggtheme = theme_classic(base_size = 11)
)$plot +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1))) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_blank()
  )

# Save SEER figures -----------------

ggsave(
  filename = "results/figures/SEER_KM_plots.pdf",
  plot = seer_km_figure,
  width = 10,
  height = 4,
  units = "in",
  device = cairo_pdf
)

############ METABRIC ###############
pred_org = readRDS("results/KM_METABRIC_19_org.rds")

pred_selected = readRDS("results/KM_METABRIC_19_selected.rds")

args1 = 19
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) *2399

p_thershold = 0.05

use_seed = 256 +start.seed
set.seed(use_seed)

trainN = paste("data/traindata","_",use_seed,".rds",sep = "")
validN = paste("data/validdata","_",use_seed,".rds",sep = "")

Ptrain = readRDS(trainN)
Pvalid = readRDS(validN)
METABRIC_data_train = importDnnetSurv(x = Ptrain[,3:71],
                                 y = Ptrain[,"OS_MONTHS"],
                                 e = Ptrain[,"OS_STATUS"])
METABRIC_data_valid = importDnnetSurv(x = Pvalid[,3:71],
                                 y = Pvalid[,"OS_MONTHS"],
                                 e = Pvalid[,"OS_STATUS"])


# Risk stratification ----------------------------------------
data_km = data.frame(y = METABRIC_data_valid@y,
                      e = METABRIC_data_valid@e,
                      strata = ifelse(pred_selected[[1]] >= median(pred_selected[[1]]),"High","Low")) %>%
  mutate(strata = factor(strata, levels = c("Low", "High")))
fit_km = survfit(Surv(y,e)~strata, data = data_km)

metabric_risk_plot <- ggsurvplot(
  fit_km,
  data = data_km,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = FALSE,
  censor = FALSE,
  size = 1,
  palette = c("#0072B2", "#D55E00"),
  legend.title = "Risk Group",
  legend.labs = c("Low", "High"),
  xlab = "Time (years)",
  ylab = "Survival probability",
  break.time.by = 60,
  ggtheme = theme_classic(base_size = 11)
)$plot +
  scale_x_continuous(labels = function(x) x / 12, breaks = seq(0, max(data_km$y, na.rm = TRUE), by = 60)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1))) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_blank()
  )
metabric_risk_plot

# NPI -------------------------------------------------------------------------

data_km <- data.frame(
  y = METABRIC_data_valid@y,
  e = METABRIC_data_valid@e,
  NPI = METABRIC_data_valid@x[, "NPI"]
) %>%
  mutate(
    strata = case_when(
      NPI < 3.4 ~ "<3.4",
      NPI >= 3.4 & NPI <= 5.4 ~ "3.4-5.4",
      NPI > 5.4 ~ ">5.4"
    ),
    strata = factor(strata, levels = c("<3.4", "3.4-5.4", ">5.4"))
  ) %>%
  filter(!is.na(strata))

fit_km <- survfit(Surv(y, e) ~ strata, data = data_km)

metabric_npi_plot <- ggsurvplot(
  fit_km,
  data = data_km,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = FALSE,
  censor = FALSE,
  size = 1,
  palette = c("#0072B2", "#D55E00", "#009E73"),
  legend.title = "NPI",
  legend.labs = c("<3.4", "3.4-5.4", ">5.4"),
  xlab = "Time (years)",
  ylab = "Survival probability",
  break.time.by = 60,
  ggtheme = theme_classic(base_size = 11)
)$plot +
  scale_x_continuous(labels = function(x) x / 12, breaks = seq(0, max(data_km$y, na.rm = TRUE), by = 60)) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 1))) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.key.width = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    plot.title = element_blank()
  )
metabric_npi_plot

# Save METABRIC figures ---------------------------------------------------------------------

metabric_km_figure <- grid.arrange(
  metabric_npi_plot,
  metabric_risk_plot,
  nrow = 1
)

ggsave(
  filename = "results/figures/METABRIC_KM_plots.pdf",
  plot = metabric_km_figure,
  width = 10,
  height = 4,
  units = "in",
  device = cairo_pdf
)
