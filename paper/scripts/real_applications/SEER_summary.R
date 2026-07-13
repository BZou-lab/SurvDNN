# SEER Cohort Characteristics and Risk Groups Summary

library(tidyverse)
library(knitr)
library(readr)

hcc_traindata <- readRDS("data/seer_traindata_unprocessed.rds")
hcc_validdata <- readRDS("data/seer_validdata_unprocessed.rds")

hcc_alldata <- rbind(hcc_traindata, hcc_validdata)

# Summary table of SEER-HCC cohort characteristics ----------------------------------------

# Variables to summarize in hcc_alldata:
# OS_STATUS (0/1 binary for event)
# OS_MONTHS (continuous follow-up time in months)
# Sex (binary variable that is either "Male" or "Female")
# Race (categorical with three levels: "White", "Black", "Asian or Pacific Islander")
# Grade (tumor grade has three ordinal levels: 1, 2, 3)
# T_Stage (ordinal with four levels: T1, T2, T3, and T4)
# N_Stage (binary variable that is either N0 or N1)
# M_Stage (binary variable that is either M0 or M1)
# time_to_trt_days (continuous time to treatment in days)
# Radiation (binary variable that is either Yes or "No/Unknown")
# Chemotherapy (binary variable that is either Yes or No)
# Fibrosis_Score (categorical variable with 3 levels: "Ishak 0-4; No to moderate fibrosis", "Ishak 5-6; Advanced/severe fibrosis", "Unknown")
# Surgery_Type (categorical variable that is either "Surgical resection" or "Tumor destruction")
# AFP_status (binary variable that is either "Positive" or "Negative")
# primary_noninvasive (binary variable that is either Yes or No)
# tumor_size (continuous variable for tumor size; units are mm)
# first_malignant_primary (binary variable that is either Yes or No)
# num_malignant (discrete count variable for the number of malignant tumors)
# Age (continuous age in years)
# Married (binary variable that is either Yes or No)
# adjusted_household_income (binary variable that is either "<$80,000" or "$80,000+")

seer_feature_order <- c(
  "Age",
  "OS_MONTHS",
  "Sex",
  "Race",
  "Grade",
  "T_Stage",
  "N_Stage",
  "M_Stage",
  "time_to_trt_days",
  "Radiation",
  "Chemotherapy",
  "Fibrosis_Score",
  "Surgery_Type",
  "AFP_status",
  "primary_noninvasive",
  "tumor_size",
  "first_malignant_primary",
  "num_malignant",
  "Married",
  "adjusted_household_income",
  "OS_STATUS"
)

seer_continuous_features <- c(
  "Age",
  "OS_MONTHS",
  "time_to_trt_days",
  "tumor_size"
)

seer_count_features <- c(
  "num_malignant"
)

seer_binary_features <- c(
  "OS_STATUS",
  "Sex",
  "N_Stage",
  "M_Stage",
  "Radiation",
  "Chemotherapy",
  "AFP_status",
  "primary_noninvasive",
  "first_malignant_primary",
  "Married",
  "adjusted_household_income"
)

seer_categorical_features <- c(
  "Race",
  "Grade",
  "T_Stage",
  "Fibrosis_Score",
  "Surgery_Type"
)

seer_feature_labels <- c(
  Age = "Age at diagnosis, years",
  OS_MONTHS = "Follow-up time, months",
  Sex = "Male sex",
  Race = "Race",
  Grade = "Tumor grade",
  T_Stage = "T stage",
  N_Stage = "N1 stage",
  M_Stage = "M1 stage",
  time_to_trt_days = "Time from diagnosis to treatment, days",
  Radiation = "Radiation received",
  Chemotherapy = "Chemotherapy received",
  Fibrosis_Score = "Fibrosis score",
  Surgery_Type = "Surgery type",
  AFP_status = "AFP positive",
  primary_noninvasive = "Non-invasive primary tumor",
  tumor_size = "Tumor size, mm",
  first_malignant_primary = "First malignant primary",
  num_malignant = "Number of malignant tumors",
  Married = "Married",
  adjusted_household_income = "Inflation-adjusted income > $80,000",
  OS_STATUS = "Overall survival event"
)

seer_binary_event_levels <- c(
  OS_STATUS = "1",
  Sex = "Male",
  N_Stage = "N1",
  M_Stage = "M1",
  Radiation = "Yes",
  Chemotherapy = "Yes",
  AFP_status = "Positive",
  primary_noninvasive = "Yes",
  first_malignant_primary = "Yes",
  Married = "Yes",
  adjusted_household_income = "$80,000+"
)

seer_categorical_level_order <- list(
  Race = c("White", "Black", "Asian or Pacific Islander"),
  Grade = c("1", "2", "3"),
  T_Stage = c("T1", "T2", "T3", "T4"),
  Fibrosis_Score = c(
    "Ishak 0-4; No to moderate fibrosis",
    "Ishak 5-6; Advanced/severe fibrosis",
    "Unknown"
  ),
  Surgery_Type = c("Tumor destruction", "Surgical resection")
)

seer_categorical_level_labels <- list(
  Grade = c(
    "1" = "I",
    "2" = "II",
    "3" = "III"
  )
)

seer_sections <- list(
  "Patient characteristics" = c(
    "Age",
    "Sex",
    "Race",
    "Married",
    "adjusted_household_income"
  ),
  "Tumor characteristics" = c(
    "Grade",
    "T_Stage",
    "N_Stage",
    "M_Stage",
    "Fibrosis_Score",
    "AFP_status",
    "primary_noninvasive",
    "tumor_size",
    "first_malignant_primary",
    "num_malignant"
  ),
  "Treatment characteristics" = c(
    "time_to_trt_days",
    "Radiation",
    "Chemotherapy",
    "Surgery_Type"
  ),
  "Outcomes" = c(
    "OS_MONTHS",
    "OS_STATUS"
  )
)

fmt_seer_mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_seer_mean_sd_integer <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sprintf("%.0f (%.0f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

fmt_seer_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.0f [%.0f, %.0f]", qs[2], qs[1], qs[3])
}

fmt_seer_number_optional_decimals <- function(x, digits = 2) {
  whole_number <- !is.na(x) & abs(x - round(x)) < sqrt(.Machine$double.eps)
  formatted <- ifelse(
    whole_number,
    sprintf("%.0f", x),
    sub("\\.?0+$", "", sprintf(paste0("%.", digits, "f"), x))
  )
  formatted[is.na(x)] <- NA_character_
  formatted
}

fmt_seer_median_iqr_optional_decimals <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf(
    "%s [%s, %s]",
    fmt_seer_number_optional_decimals(qs[2], digits),
    fmt_seer_number_optional_decimals(qs[1], digits),
    fmt_seer_number_optional_decimals(qs[3], digits)
  )
}

fmt_seer_n_pct <- function(x, level) {
  x_chr <- as.character(x)
  n_nonmissing <- sum(!is.na(x_chr))
  n_level <- sum(x_chr == level, na.rm = TRUE)
  if (n_nonmissing == 0) {
    return("0 (NA)")
  }
  if (n_level < 11) {
    return("<11")
  }
  sprintf("%d (%.1f%%)", n_level, 100 * n_level / n_nonmissing)
}

get_seer_analysis_variable <- function(data, feature) {
  x <- data[[feature]]
  x
}

make_seer_summary_rows <- function(data, feature) {
  if (feature %in% seer_continuous_features) {
    x <- get_seer_analysis_variable(data, feature)
    summary_value <- if (feature == "OS_MONTHS") {
      fmt_seer_median_iqr_optional_decimals(x)
    } else if (feature %in% c("time_to_trt_days", "tumor_size")) {
      fmt_seer_median_iqr(x)
    } else {
      fmt_seer_median_iqr_optional_decimals(x)
    }
    
    return(tibble(
      Feature = feature,
      Characteristic = unname(seer_feature_labels[feature]),
      Overall = summary_value
    ))
  }
  
  if (feature %in% seer_count_features) {
    return(tibble(
      Feature = feature,
      Characteristic = unname(seer_feature_labels[feature]),
      Overall = fmt_seer_median_iqr(data[[feature]])
    ))
  }
  
  if (feature %in% seer_binary_features) {
    return(tibble(
      Feature = feature,
      Characteristic = unname(seer_feature_labels[feature]),
      Overall = fmt_seer_n_pct(data[[feature]], unname(seer_binary_event_levels[feature]))
    ))
  }
  
  if (feature %in% seer_categorical_features) {
    observed_levels <- unique(as.character(data[[feature]][!is.na(data[[feature]])]))
    ordered_levels <- seer_categorical_level_order[[feature]]
    ordered_levels <- c(ordered_levels[ordered_levels %in% observed_levels],
                        setdiff(sort(observed_levels), ordered_levels))
    display_levels <- ordered_levels
    if (!is.null(seer_categorical_level_labels[[feature]])) {
      mapped_levels <- unname(seer_categorical_level_labels[[feature]][ordered_levels])
      display_levels <- ifelse(is.na(mapped_levels), ordered_levels, mapped_levels)
    }
    
    return(bind_rows(
      tibble(
        Feature = feature,
        Characteristic = unname(seer_feature_labels[feature]),
        Overall = ""
      ),
      tibble(
        Feature = feature,
        Characteristic = paste0("  ", display_levels),
        Overall = map_chr(ordered_levels, ~ fmt_seer_n_pct(data[[feature]], .x))
      )
    ))
  }
  
  tibble(
    Feature = feature,
    Characteristic = unname(seer_feature_labels[feature]),
    Overall = NA_character_
  )
}

make_seer_summary_table <- function(hcc_alldata) {
  map_dfr(seer_feature_order, ~ make_seer_summary_rows(hcc_alldata, .x))
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

make_seer_summary_latex <- function(hcc_alldata) {
  summary_table <- make_seer_summary_table(hcc_alldata)
  n_label <- paste0("Overall (N = ", format(nrow(hcc_alldata), big.mark = ","), ")")
  
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\caption{SEER-HCC cohort characteristics}",
    "\\begin{tabular}{ll}",
    "\\toprule",
    paste0("Characteristic & ", latex_escape(n_label), " \\\\"),
    "\\midrule"
  )
  
  for (section_name in names(seer_sections)) {
    section_rows <- summary_table %>%
      filter(Feature %in% seer_sections[[section_name]]) %>%
      mutate(Feature = factor(Feature, levels = seer_sections[[section_name]])) %>%
      arrange(Feature) %>%
      dplyr::select(-Feature)
    
    lines <- c(
      lines,
      paste0("\\multicolumn{2}{l}{\\textbf{", latex_escape(section_name), "}} \\\\")
    )
    
    for (i in seq_len(nrow(section_rows))) {
      characteristic <- section_rows$Characteristic[i]
      if (startsWith(characteristic, "  ")) {
        characteristic <- paste0("\\hspace{1em}", trimws(characteristic))
      } else {
        characteristic <- latex_escape(characteristic)
      }
      
      lines <- c(
        lines,
        paste0(characteristic, " & ", latex_escape(section_rows$Overall[i]), " \\\\")
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

seer_summary_table <- make_seer_summary_table(hcc_alldata)
seer_summary_latex <- make_seer_summary_latex(hcc_alldata)

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
writeLines(seer_summary_latex, "results/tables/seer_hcc_cohort_characteristics.tex")


# Summary table of SEER risk group characteristics ----------------------------------------------

# Stratify subjects into risk groups 
risk_predictions = readRDS("results/KM_SEER_71_selected.rds")
strata <- ifelse(risk_predictions[[1]] >= median(risk_predictions[[1]], na.rm = TRUE), "High", "Low")

# Add row for risk group to the test set 
hcc_validdata$risk_group <- strata[rownames(hcc_validdata)]
if (all(is.na(hcc_validdata$risk_group)) && length(strata) == nrow(hcc_validdata)) {
  hcc_validdata$risk_group <- strata
}

# Create risk group comparison table using the same variables as the cohort summary table.
seer_risk_group_feature_order <- seer_feature_order

seer_risk_group_continuous_vars <- seer_continuous_features

seer_risk_group_count_vars <- seer_count_features

seer_risk_group_binary_vars <- seer_binary_features

seer_risk_group_categorical_vars <- seer_categorical_features

format_seer_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return(NA_character_)
  }
  if (p_value < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p_value)
}

test_seer_continuous_by_risk_group <- function(data, variable) {
  x <- get_seer_analysis_variable(data, variable)
  test_data <- tibble(
    value = suppressWarnings(as.numeric(x)),
    risk_group = data$risk_group
  ) %>%
    filter(!is.na(value), !is.na(risk_group))
  
  if (n_distinct(test_data$risk_group) < 2) {
    return(NA_real_)
  }
  
  wilcox.test(value ~ risk_group, data = test_data)$p.value
}

test_seer_categorical_by_risk_group <- function(data, variable) {
  test_data <- data %>%
    filter(!is.na(.data[[variable]]), !is.na(risk_group))
  tab <- table(test_data[[variable]], test_data$risk_group)
  
  if (any(dim(tab) < 2)) {
    return(NA_real_)
  }
  
  chi_result <- suppressWarnings(chisq.test(tab, correct = FALSE))
  if (any(chi_result$expected < 5)) {
    tryCatch(
      fisher.test(tab)$p.value,
      error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value
    )
  } else {
    chi_result$p.value
  }
}

format_seer_median_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.2f [%.2f, %.2f]", qs[2], qs[1], qs[3])
}

format_seer_median_iqr_integer <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  sprintf("%.0f [%.0f, %.0f]", qs[2], qs[1], qs[3])
}

format_seer_group_n_pct <- function(x, level) {
  fmt_seer_n_pct(x, level)
}

make_seer_risk_group_continuous_row <- function(data, variable) {
  low_data <- data %>% filter(risk_group == "Low")
  high_data <- data %>% filter(risk_group == "High")
  formatter <- if (variable == "OS_MONTHS") {
    fmt_seer_median_iqr_optional_decimals
  } else if (variable %in% c("time_to_trt_days", "tumor_size", "num_malignant")) {
    format_seer_median_iqr_integer
  } else {
    format_seer_median_iqr
  }
  
  tibble(
    Feature = variable,
    Characteristic = unname(seer_feature_labels[variable]),
    `Low Risk` = formatter(get_seer_analysis_variable(low_data, variable)),
    `High Risk` = formatter(get_seer_analysis_variable(high_data, variable)),
    `P-value` = format_seer_p_value(test_seer_continuous_by_risk_group(data, variable))
  )
}

make_seer_risk_group_binary_row <- function(data, variable) {
  low_data <- data %>% filter(risk_group == "Low")
  high_data <- data %>% filter(risk_group == "High")
  event_level <- unname(seer_binary_event_levels[variable])
  
  tibble(
    Feature = variable,
    Characteristic = unname(seer_feature_labels[variable]),
    `Low Risk` = format_seer_group_n_pct(low_data[[variable]], event_level),
    `High Risk` = format_seer_group_n_pct(high_data[[variable]], event_level),
    `P-value` = format_seer_p_value(test_seer_categorical_by_risk_group(data, variable))
  )
}

make_seer_risk_group_categorical_rows <- function(data, variable) {
  low_data <- data %>% filter(risk_group == "Low")
  high_data <- data %>% filter(risk_group == "High")
  observed_levels <- unique(as.character(data[[variable]][!is.na(data[[variable]])]))
  ordered_levels <- seer_categorical_level_order[[variable]]
  ordered_levels <- c(ordered_levels[ordered_levels %in% observed_levels],
                      setdiff(sort(observed_levels), ordered_levels))
  display_levels <- ordered_levels
  if (!is.null(seer_categorical_level_labels[[variable]])) {
    mapped_levels <- unname(seer_categorical_level_labels[[variable]][ordered_levels])
    display_levels <- ifelse(is.na(mapped_levels), ordered_levels, mapped_levels)
  }
  
  bind_rows(
    tibble(
      Feature = variable,
      Characteristic = unname(seer_feature_labels[variable]),
      `Low Risk` = "",
      `High Risk` = "",
      `P-value` = format_seer_p_value(test_seer_categorical_by_risk_group(data, variable))
    ),
    tibble(
      Feature = variable,
      Characteristic = paste0("  ", display_levels),
      `Low Risk` = map_chr(ordered_levels, ~ format_seer_group_n_pct(low_data[[variable]], .x)),
      `High Risk` = map_chr(ordered_levels, ~ format_seer_group_n_pct(high_data[[variable]], .x)),
      `P-value` = ""
    )
  )
}

make_seer_risk_group_rows <- function(data, variable) {
  if (variable %in% seer_risk_group_continuous_vars) {
    return(make_seer_risk_group_continuous_row(data, variable))
  }
  if (variable %in% seer_risk_group_count_vars) {
    return(make_seer_risk_group_continuous_row(data, variable))
  }
  if (variable %in% seer_risk_group_binary_vars) {
    return(make_seer_risk_group_binary_row(data, variable))
  }
  if (variable %in% seer_risk_group_categorical_vars) {
    return(make_seer_risk_group_categorical_rows(data, variable))
  }
  
  tibble(
    Feature = variable,
    Characteristic = unname(seer_feature_labels[variable]),
    `Low Risk` = NA_character_,
    `High Risk` = NA_character_,
    `P-value` = NA_character_
  )
}

seer_risk_group_data <- hcc_validdata %>%
  mutate(risk_group = factor(risk_group, levels = c("Low", "High")))

seer_risk_group_summary_df <- map_dfr(
  seer_risk_group_feature_order,
  ~ make_seer_risk_group_rows(seer_risk_group_data, .x)
) %>%
  rename(
    !!paste0("Low Risk (N = ", sum(seer_risk_group_data$risk_group == "Low", na.rm = TRUE), ")") := `Low Risk`,
    !!paste0("High Risk (N = ", sum(seer_risk_group_data$risk_group == "High", na.rm = TRUE), ")") := `High Risk`
  )

make_seer_risk_group_latex <- function(summary_table) {
  low_col <- grep("^Low Risk", names(summary_table), value = TRUE)
  high_col <- grep("^High Risk", names(summary_table), value = TRUE)
  
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\small",
    "\\caption{SEER-HCC validation cohort characteristics by SurvDNN risk group}",
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
  
  for (section_name in names(seer_sections)) {
    section_rows <- summary_table %>%
      filter(Feature %in% seer_sections[[section_name]]) %>%
      mutate(Feature = factor(Feature, levels = seer_sections[[section_name]])) %>%
      arrange(Feature) %>%
      dplyr::select(-Feature)
    
    lines <- c(
      lines,
      paste0("\\multicolumn{4}{l}{\\textbf{", latex_escape(section_name), "}} \\\\")
    )
    
    for (i in seq_len(nrow(section_rows))) {
      characteristic <- section_rows$Characteristic[i]
      if (startsWith(characteristic, "  ")) {
        characteristic <- paste0("\\hspace{1em}", trimws(characteristic))
      } else {
        characteristic <- latex_escape(characteristic)
      }
      
      lines <- c(
        lines,
        paste0(
          characteristic,
          " & ",
          latex_escape(section_rows[[low_col]][i]),
          " & ",
          latex_escape(section_rows[[high_col]][i]),
          " & ",
          latex_escape(section_rows$`P-value`[i]),
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

seer_risk_group_latex <- make_seer_risk_group_latex(seer_risk_group_summary_df)
writeLines(seer_risk_group_latex, "results/tables/seer_hcc_risk_group_characteristics.tex")

















