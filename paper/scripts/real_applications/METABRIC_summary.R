# METABRIC Cohort Characteristics and Risk Groups Summary

library(tidyverse)
library(knitr)

metabric_traindata <- readRDS("data/traindata_45837.rds")
metabric_validdata <- readRDS("data/validdata_45837.rds")

metabric_alldata <- rbind(metabric_traindata, metabric_validdata)

# Summary table of METABRIC cohort characteristics ----------------------------------------

# Variables to summarize in metabric_alldata:
# OS_STATUS (0/1 binary for event)
# OS_MONTHS (continuous follow-up time, converted from months to years)
# AGE_AT_DIAGNOSIS (continuous)
# LYMPH_NODES_EXAMINED_POSITIVE (discrete count)
# RFS_STATUS_R (0/1 binary)
# ER_STATUS_P (0/1 binary)
# NPI (continuous)
# TP53 (0/1 binary)
# PR_STATUS_P (0/1 binary)
# TUMOR_SIZE (continuous)
# CHEMOTHERAPY_Y (0/1 binary)
# HORMONE_THERAPY_Y (0/1 binary)
# HER2_STATUS_P (0/1 binary)
# GATA3 (0/1 binary)
# NOTCH1 (0/1 binary)
# SETD1A (0/1 binary)
# EP300 (0/1 binary)
# COL12A1 (0/1 binary)
# THADA (0/1 binary)
# CBFB (0/1 binary)


feature_order <- c(
  "AGE_AT_DIAGNOSIS",
  "LYMPH_NODES_EXAMINED_POSITIVE",
  "NPI",
  "TUMOR_SIZE",
  "ER_STATUS_P",
  "PR_STATUS_P",
  "HER2_STATUS_P",
  "TP53",
  "GATA3",
  "NOTCH1",
  "SETD1A",
  "EP300",
  "COL12A1",
  "THADA",
  "CBFB",
  "CHEMOTHERAPY_Y",
  "HORMONE_THERAPY_Y",
  "OS_MONTHS",
  "OS_STATUS",
  "RFS_STATUS_R"
)

continuous_features <- c(
  "OS_MONTHS",
  "AGE_AT_DIAGNOSIS",
  "NPI",
  "TUMOR_SIZE"
)

count_features <- c(
  "LYMPH_NODES_EXAMINED_POSITIVE"
)

binary_features <- c(
  "OS_STATUS",
  "RFS_STATUS_R",
  "ER_STATUS_P",
  "TP53",
  "GATA3",
  "NOTCH1",
  "SETD1A",
  "EP300",
  "COL12A1",
  "THADA",
  "CBFB",
  "PR_STATUS_P",
  "CHEMOTHERAPY_Y",
  "HORMONE_THERAPY_Y",
  "HER2_STATUS_P"
)

feature_labels <- c(
  OS_STATUS = "Overall survival event",
  OS_MONTHS = "Follow-up time, years",
  AGE_AT_DIAGNOSIS = "Age at diagnosis, years",
  LYMPH_NODES_EXAMINED_POSITIVE = "Positive lymph nodes examined",
  RFS_STATUS_R = "Recurrence-free survival event",
  ER_STATUS_P = "ER positive",
  NPI = "Nottingham Prognostic Index",
  TP53 = "TP53 mutation",
  GATA3 = "GATA3 mutation",
  NOTCH1 = "NOTCH1 mutation",
  SETD1A = "SETD1A mutation",
  EP300 = "EP300 mutation",
  COL12A1 = "COL12A1 mutation",
  THADA = "THADA mutation",
  CBFB = "CBFB mutation",
  PR_STATUS_P = "PR positive",
  TUMOR_SIZE = "Tumor size, mm",
  CHEMOTHERAPY_Y = "Chemotherapy received",
  HORMONE_THERAPY_Y = "Hormone therapy received",
  HER2_STATUS_P = "HER2 positive"
)

fmt_mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_mean_sd_integer <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sprintf("%.0f (%.0f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.0f [%.0f, %.0f]", qs[2], qs[1], qs[3])
}

fmt_number_optional_decimals <- function(x, digits = 2) {
  whole_number <- !is.na(x) & abs(x - round(x)) < sqrt(.Machine$double.eps)
  formatted <- ifelse(
    whole_number,
    sprintf("%.0f", x),
    sprintf(paste0("%.", digits, "f"), x)
  )
  formatted[is.na(x)] <- NA_character_
  formatted
}

fmt_median_iqr_optional_decimals <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(
    "%s [%s, %s]",
    fmt_number_optional_decimals(qs[2], digits),
    fmt_number_optional_decimals(qs[1], digits),
    fmt_number_optional_decimals(qs[3], digits)
  )
}

fmt_binary <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  n_nonmissing <- sum(!is.na(x))
  n_yes <- sum(x == 1, na.rm = TRUE)
  if (n_nonmissing == 0) {
    return("0 (NA)")
  }
  if (n_yes < 11) {
    return("<11")
  }
  pct_yes <- 100 * n_yes / n_nonmissing
  sprintf("%d (%.1f%%)", n_yes, pct_yes)
}

make_summary_rows <- function(data, features, formatter) {
  tibble(
    Characteristic = unname(feature_labels[features]),
    Overall = map_chr(features, function(feature) {
      x <- data[[feature]]
      if (feature == "OS_MONTHS") {
        x <- x / 12
      }
      formatter(x)
    })
  )
}

make_metabric_summary_table <- function(metabric_alldata) {
  bind_rows(
    make_summary_rows(metabric_alldata, continuous_features, fmt_median_iqr_optional_decimals),
    make_summary_rows(metabric_alldata, count_features, fmt_median_iqr),
    make_summary_rows(metabric_alldata, binary_features, fmt_binary)
  ) %>%
    mutate(
      Characteristic = factor(
        Characteristic,
        levels = unname(feature_labels[feature_order])
      )
    ) %>%
    arrange(Characteristic) %>%
    mutate(Characteristic = as.character(Characteristic))
}

make_metabric_summary_kable <- function(metabric_alldata) {
  summary_table <- make_metabric_summary_table(metabric_alldata)
  
  kable(
    summary_table,
    format = "latex",
    caption = "METABRIC cohort characteristics",
    col.names = c("Characteristic", paste0("Overall (N = ", format(nrow(metabric_alldata), big.mark = ","), ")")),
    align = c("l", "l")
  )
}

latex_escape <- function(x) {
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("#", "\\#", x, fixed = TRUE)
  x <- gsub("<", "\\textless{}", x, fixed = TRUE)
  x
}

make_metabric_summary_latex <- function(metabric_alldata) {
  summary_table <- make_metabric_summary_table(metabric_alldata)
  n_label <- paste0("Overall (N = ", format(nrow(metabric_alldata), big.mark = ","), ")")
  
  section_rows <- list(
    "Patient characteristics" = c(
      "Age at diagnosis, years"
    ),
    "Tumor characteristics" = c(
      "Positive lymph nodes examined",
      "Nottingham Prognostic Index",
      "Tumor size, mm",
      "ER positive",
      "PR positive",
      "HER2 positive",
      "TP53 mutation",
      "GATA3 mutation",
      "NOTCH1 mutation",
      "SETD1A mutation",
      "EP300 mutation",
      "COL12A1 mutation",
      "THADA mutation",
      "CBFB mutation"
    ),
    "Treatment characteristics" = c(
      "Chemotherapy received",
      "Hormone therapy received"
    ),
    "Outcomes" = c(
      "Follow-up time, years",
      "Overall survival event",
      "Recurrence-free survival event"
    )
  )
  
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{METABRIC cohort characteristics}",
    "\\begin{tabular}{ll}",
    "\\toprule",
    paste0("Characteristic & ", latex_escape(n_label), " \\\\"),
    "\\midrule"
  )
  
  for (section_name in names(section_rows)) {
    rows <- summary_table %>%
      filter(Characteristic %in% section_rows[[section_name]]) %>%
      mutate(
        Characteristic = factor(Characteristic, levels = section_rows[[section_name]])
      ) %>%
      arrange(Characteristic) %>%
      mutate(Characteristic = as.character(Characteristic))
    
    lines <- c(
      lines,
      paste0("\\multicolumn{2}{l}{\\textbf{", latex_escape(section_name), "}} \\\\")
    )
    
    for (i in seq_len(nrow(rows))) {
      lines <- c(
        lines,
        paste0(latex_escape(rows$Characteristic[i]), " & ", latex_escape(rows$Overall[i]), " \\\\")
      )
    }
    
    lines <- c(lines, "\\midrule")
  }
  
  lines[length(lines)] <- "\\bottomrule"
  
  c(
    lines,
    "\\end{tabular}",
    "\\begin{minipage}{0.95\\linewidth}",
    "\\vspace{0.5em}",
    "\\footnotesize Continuous variables are presented as median [IQR]. Categorical variables are presented as n (\\%); cells with counts below 11 are shown as <11 without percentages.",
    "\\end{minipage}",
    "\\end{table}"
  )
}

metabric_summary_latex <- make_metabric_summary_latex(metabric_alldata)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
writeLines(metabric_summary_latex, "results/tables/metabric_cohort_characteristics.tex")


# Summary table of METABRIC risk group characteristics ----------------------------------------------

# Stratify subjects into risk groups 
risk_predictions = readRDS("results/KM_METABRIC_19_selected.rds")
strata <- ifelse(risk_predictions[[1]] >= median(risk_predictions[[1]]),"High","Low")

# Add row for risk group to the test set 
metabric_validdata$risk_group <- strata[rownames(metabric_validdata)]

# Create risk group comparison table 
risk_group_continuous_vars <- c(
  "AGE_AT_DIAGNOSIS",
  "FOLLOW_UP_YEARS",
  "LYMPH_NODES_EXAMINED_POSITIVE",
  "NPI",
  "TUMOR_SIZE"
)

risk_group_binary_vars <- c(
  "ER_STATUS_P",
  "PR_STATUS_P",
  "HER2_STATUS_P",
  "TP53",
  "GATA3",
  "NOTCH1",
  "SETD1A",
  "EP300",
  "COL12A1",
  "THADA",
  "CBFB",
  "CHEMOTHERAPY_Y",
  "HORMONE_THERAPY_Y",
  "OS_STATUS",
  "RFS_STATUS_R"
)

risk_group_feature_order <- c(
  "AGE_AT_DIAGNOSIS",
  "LYMPH_NODES_EXAMINED_POSITIVE",
  "NPI",
  "TUMOR_SIZE",
  "ER_STATUS_P",
  "PR_STATUS_P",
  "HER2_STATUS_P",
  "TP53",
  "GATA3",
  "NOTCH1",
  "SETD1A",
  "EP300",
  "COL12A1",
  "THADA",
  "CBFB",
  "CHEMOTHERAPY_Y",
  "HORMONE_THERAPY_Y",
  "FOLLOW_UP_YEARS",
  "OS_STATUS",
  "RFS_STATUS_R"
)

risk_group_feature_labels <- c(
  AGE_AT_DIAGNOSIS = "Age at diagnosis, years",
  FOLLOW_UP_YEARS = "Follow-up time, years",
  LYMPH_NODES_EXAMINED_POSITIVE = "Positive lymph nodes examined",
  NPI = "Nottingham Prognostic Index",
  TUMOR_SIZE = "Tumor size, mm",
  ER_STATUS_P = "ER positive",
  PR_STATUS_P = "PR positive",
  HER2_STATUS_P = "HER2 positive",
  TP53 = "TP53 mutation",
  GATA3 = "GATA3 mutation",
  NOTCH1 = "NOTCH1 mutation",
  SETD1A = "SETD1A mutation",
  EP300 = "EP300 mutation",
  COL12A1 = "COL12A1 mutation",
  THADA = "THADA mutation",
  CBFB = "CBFB mutation",
  CHEMOTHERAPY_Y = "Chemotherapy received",
  HORMONE_THERAPY_Y = "Hormone therapy received",
  OS_STATUS = "Overall survival event",
  RFS_STATUS_R = "Recurrence-free survival event"
)

risk_group_sections <- list(
  "Patient characteristics" = c("AGE_AT_DIAGNOSIS"),
  "Tumor characteristics" = c(
    "LYMPH_NODES_EXAMINED_POSITIVE",
    "NPI",
    "TUMOR_SIZE",
    "ER_STATUS_P",
    "PR_STATUS_P",
    "HER2_STATUS_P",
    "TP53",
    "GATA3",
    "NOTCH1",
    "SETD1A",
    "EP300",
    "COL12A1",
    "THADA",
    "CBFB"
  ),
  "Treatment characteristics" = c("CHEMOTHERAPY_Y", "HORMONE_THERAPY_Y"),
  "Outcomes" = c("FOLLOW_UP_YEARS", "OS_STATUS", "RFS_STATUS_R")
)

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return(NA_character_)
  }
  if (p_value < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p_value)
}

test_binary_by_risk_group <- function(data, variable) {
  tab <- table(data[[variable]], data$risk_group)
  
  if (any(dim(tab) < 2)) {
    return(NA_real_)
  }
  
  chi_result <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (any(chi_result$expected < 5)) {
    fisher.test(tab)$p.value
  } else {
    chi_result$p.value
  }
}

test_continuous_by_risk_group <- function(data, variable) {
  wilcox.test(data[[variable]] ~ data$risk_group)$p.value
}

format_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.2f [%.2f, %.2f]", qs[2], qs[1], qs[3])
}

format_median_iqr_integer <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.0f [%.0f, %.0f]", qs[2], qs[1], qs[3])
}

format_yes_n_pct <- function(x) {
  n_nonmissing <- sum(!is.na(x))
  n_yes <- sum(x == "Yes" | x == 1, na.rm = TRUE)
  if (n_nonmissing == 0) {
    return("0 (NA)")
  }
  if (n_yes < 11) {
    return("<11")
  }
  sprintf("%d (%.1f%%)", n_yes, 100 * n_yes / n_nonmissing)
}

make_risk_group_manual_row <- function(data, variable) {
  low_data <- data %>% filter(risk_group == "Low")
  high_data <- data %>% filter(risk_group == "High")
  
  if (variable %in% risk_group_continuous_vars) {
    formatter <- if (variable == "LYMPH_NODES_EXAMINED_POSITIVE") {
      format_median_iqr_integer
    } else {
      fmt_median_iqr_optional_decimals
    }
    low_summary <- formatter(low_data[[variable]])
    high_summary <- formatter(high_data[[variable]])
    p_value <- test_continuous_by_risk_group(data, variable)
  } else {
    low_summary <- format_yes_n_pct(low_data[[variable]])
    high_summary <- format_yes_n_pct(high_data[[variable]])
    p_value <- test_binary_by_risk_group(data, variable)
  }
  
  tibble(
    Characteristic = unname(risk_group_feature_labels[variable]),
    `Low Risk` = low_summary,
    `High Risk` = high_summary,
    `P-value` = format_p_value(p_value)
  )
}

metabric_risk_group_data <- metabric_validdata %>%
  mutate(
    risk_group = factor(risk_group, levels = c("Low", "High")),
    FOLLOW_UP_YEARS = OS_MONTHS / 12,
    across(all_of(risk_group_binary_vars), ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes")))
  )

metabric_risk_group_summary_df <- map_dfr(
  risk_group_feature_order,
  ~ make_risk_group_manual_row(metabric_risk_group_data, .x)
) %>%
  rename(
    !!paste0("Low Risk (N = ", sum(metabric_risk_group_data$risk_group == "Low", na.rm = TRUE), ")") := `Low Risk`,
    !!paste0("High Risk (N = ", sum(metabric_risk_group_data$risk_group == "High", na.rm = TRUE), ")") := `High Risk`
  )

# Create a booktabs-style LaTeX table matching the cohort summary table.
make_metabric_risk_group_latex <- function(summary_table) {
  low_col <- grep("^Low Risk", names(summary_table), value = TRUE)
  high_col <- grep("^High Risk", names(summary_table), value = TRUE)
  
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{METABRIC validation cohort characteristics by SurvDNN risk group}",
    "\\begin{tabular}{llll}",
    "\\toprule",
    paste0(
      "Characteristic & ",
      latex_escape(low_col),
      " & ",
      latex_escape(high_col),
      " & P-value \\\\"
    ),
    "\\midrule"
  )
  
  for (section_name in names(risk_group_sections)) {
    section_labels <- unname(risk_group_feature_labels[risk_group_sections[[section_name]]])
    rows <- summary_table %>%
      filter(Characteristic %in% section_labels) %>%
      mutate(Characteristic = factor(Characteristic, levels = section_labels)) %>%
      arrange(Characteristic) %>%
      mutate(Characteristic = as.character(Characteristic))
    
    lines <- c(
      lines,
      paste0("\\multicolumn{4}{l}{\\textbf{", latex_escape(section_name), "}} \\\\")
    )
    
    for (i in seq_len(nrow(rows))) {
      lines <- c(
        lines,
        paste0(
          latex_escape(rows$Characteristic[i]),
          " & ",
          latex_escape(rows[[low_col]][i]),
          " & ",
          latex_escape(rows[[high_col]][i]),
          " & ",
          latex_escape(rows$`P-value`[i]),
          " \\\\"
        )
      )
    }
    
    lines <- c(lines, "\\midrule")
  }
  
  lines[length(lines)] <- "\\bottomrule"
  
  c(
    lines,
    "\\end{tabular}",
    "\\begin{minipage}{0.95\\linewidth}",
    "\\vspace{0.5em}",
    "\\footnotesize Continuous variables are presented as median [IQR]. Categorical variables are presented as n (\\%); cells with counts below 11 are shown as <11 without percentages. P-values are from Wilcoxon rank-sum tests for continuous variables and chi-square or Fisher's exact tests for categorical variables.",
    "\\end{minipage}",
    "\\end{table}"
  )
}

metabric_risk_group_latex <- make_metabric_risk_group_latex(metabric_risk_group_summary_df)
writeLines(metabric_risk_group_latex, "results/tables/metabric_risk_group_characteristics.tex")







