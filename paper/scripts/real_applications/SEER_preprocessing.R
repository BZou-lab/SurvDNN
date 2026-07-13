# Script to clean and preprocess SEER HCC data for survival analysis

rm(list = ls())
library(tidyverse)
library(data.table)

df = fread("data/HCC1.txt") %>% as.data.frame()

# Extract cohort and code variables according to SEER data dictionary -------------------------

df_sub = df %>% filter(V4 == "2020" | V4 == "2021" | V4 == "2019" | V4 == "2018") %>%
  filter(V5 != "American Indian/Alaska Native" & V5 != "Unknown") %>% #Race
  filter(V6 == 1 | V6 == 2|V6 == 3) %>% # Derived Grade 
  filter(V11 == "T1" | V11 == "T1a" | V11 == "T1b" | V11 == "T2" | V11 == "T3" |V11 == "T4" ) %>% # T stage
  filter(V12 == "N0" | V12 == "N1") %>% # N stage
  filter(V13 == "M0" | V13 == "M1") %>% # M stage 
  filter(V14 != "731+ days" & V14 !="Unable to calculate") %>% # time trt
  filter(V15 != "Recommended, unknown if administered" & V15 != "Refused (1988+)") %>% # Radiation
  filter(V18 != 90 & V18 != 99) %>%  # 0;10-19;20-80 
  filter(V20 == "Negative/normal; within normal limits" | V20 == "Positive/elevated") %>% # AFP
  filter(V27 != "999" & V27 != "998" & V27 != "990") %>% # Tumor Size
  filter(V2 != "01-04 years" & V2 != "05-09 years" & V2 != "10-14 years" & V2 != "15-19 years") %>%
  filter(V38 != "Unknown")

df_sub = df_sub[,c("V1","V3","V5","V6","V11","V12","V13","V14","V15","V16","V18","V19","V20","V24","V27","V32","V33","V34","V35","V37","V38","V39")]

df_sub$Sex_F = ifelse(df_sub$V3 == "Female",1,0)

df_sub$race_black = ifelse(df_sub$V5 == "Black",1,0)
df_sub$race_API = ifelse(df_sub$V5 =="Asian or Pacific Islander",1,0)

df_sub$grade_2 = ifelse(df_sub$V6 == 2,1,0)
df_sub$grade_3 = ifelse(df_sub$V6 == 3,1,0)  

df_sub$T2 = ifelse(df_sub$V11 =="T2",1,0) 
df_sub$T3 = ifelse(df_sub$V11 =="T3",1,0)
df_sub$T4 = ifelse(df_sub$V11 =="T4",1,0)

df_sub$N1 = ifelse(df_sub$V12 =="N1",1,0) 

df_sub$M1 = ifelse(df_sub$V13 =="M1",1,0)

df_sub$diag_to_trt = as.numeric(df_sub$V14)

df_sub$Radiation_Y = ifelse(df_sub$V15 != "None/Unknown",1,0)

df_sub$Chemotherapy_Y = ifelse(df_sub$V16 == "Yes",1,0)

df_sub$V18 = as.numeric(df_sub$V18)

df_sub$surg_destruction = ifelse(df_sub$V18 >= 10 & df_sub$V18 < 20,1,0)
df_sub$surg_resection = ifelse(df_sub$V18 >= 21 & df_sub$V18 < 90,1,0)

df_sub$fibrosis_low = ifelse(df_sub$V19 == "Ishak 0-4; No to moderate fibrosis; METAVIR F0-F3; Batt-Ludwig 0-3",1,0)
df_sub$fibrosis_high = ifelse(df_sub$V19 == "Ishak 5-6; Advanced/severe fibrosis; METAVIR F4; Batt-Ludwig 4; Cirrhosis",1,0)

df_sub$AFP_pos = ifelse(df_sub$V20 == "Positive/elevated",1,0)

df_sub$primary_noninvasive_Y = ifelse(df_sub$V24 == "100",1,0)

df_sub$tumor_size = as.numeric(df_sub$V27)

df_sub$OS_STATUS = ifelse(df_sub$V32 == "Dead",1,0)
df_sub$OS_MONTHS = as.numeric(df_sub$V33)

df_sub$first_malignant_primary = ifelse(df_sub$V34 == "Yes",1,0)

df_sub$num_malignant = as.numeric(df_sub$V35)

df_sub$age = apply(df_sub,1,function(x){
  age_str = x["V37"]
  age_str2 = str_split(age_str," ")[[1]][1]
  if (age_str2 == "90+"){
    return(90)
  }else{
    age = as.numeric(str_split(age_str," ")[[1]][1])
    return(age)
  }
})

df_sub$Married = ifelse(df_sub$V38 == "Married (including common law)",1,0)

df_sub$Income_high = ifelse(df_sub$V39 %in% c("< $40,000", "$40,000 - $44,999", "$45,000 - $49,999", 
                                      "$50,000 - $54,999", "$55,000 - $59,999", "$60,000 - $64,999",
                                      "$65,000 - $69,999", "$70,000 - $74,999", "$75,000 - $79,999"),0,1)


df_sub = df_sub[,23:49]
library(survival)

summary(coxph(Surv(OS_MONTHS,OS_STATUS)~.,data = df_sub))

pheatmap::pheatmap(cor(df_sub))
  
# Save final cohort dataset for modeling purposes 
saveRDS(df_sub,file = "data/SEER_HCC_2018_2021.rds")



# Save copy of raw cohort data before one-hot encoding ---------------------------------------------
df_raw = df %>% filter(V4 == "2020" | V4 == "2021" | V4 == "2019" | V4 == "2018") %>%
  filter(V5 != "American Indian/Alaska Native" & V5 != "Unknown") %>% #Race
  filter(V6 == 1 | V6 == 2|V6 == 3) %>% # Derived Grade 
  filter(V11 == "T1" | V11 == "T1a" | V11 == "T1b" | V11 == "T2" | V11 == "T3" |V11 == "T4" ) %>% # T stage
  filter(V12 == "N0" | V12 == "N1") %>% # N stage
  filter(V13 == "M0" | V13 == "M1") %>% # M stage 
  filter(V14 != "731+ days" & V14 !="Unable to calculate") %>% # time trt
  filter(V15 != "Recommended, unknown if administered" & V15 != "Refused (1988+)") %>% # Radiation
  filter(V18 != 90 & V18 != 99) %>%  # 0;10-19;20-80 
  filter(V20 == "Negative/normal; within normal limits" | V20 == "Positive/elevated") %>% # AFP
  filter(V27 != "999" & V27 != "998" & V27 != "990") %>% # Tumor Size
  filter(V2 != "01-04 years" & V2 != "05-09 years" & V2 != "10-14 years" & V2 != "15-19 years") %>%
  filter(V38 != "Unknown")


df_raw = df_raw[,c("V3","V5","V6","V11","V12","V13","V14","V15","V16","V18","V19","V20","V24","V27","V32","V33","V34","V35","V37","V38","V39")]

df_raw$Sex <- df_raw$V3
df_raw$Race <- df_raw$V5
df_raw$Grade <- df_raw$V6
df_raw <- df_raw %>% mutate(T_Stage = case_when(V11 %in% c("T1a", "T1b") ~ "T1",
                                                V11 == "T2" ~ "T2",
                                                V11 == "T3" ~ "T3",
                                                V11 == "T4" ~ "T4"))
df_raw$N_Stage <- df_raw$V12
df_raw$M_Stage <- df_raw$V13
df_raw$time_to_trt_days <- as.numeric(df_raw$V14)
df_raw <- df_raw %>% mutate(Radiation = case_when(V15 == "None/Unknown" ~ "No/Unknown",
                                                            TRUE ~ "Yes"))

df_raw <- df_raw %>% mutate(Chemotherapy = case_when(V16 == "No/Unknown" ~ "No/Unknown",
                                                               TRUE ~ "Yes"))

df_raw <- df_raw %>% mutate(Surgery_Type = case_when((V18 >= 10 & V18 < 20) ~ "Tumor destruction",
                            TRUE ~ "Surgical resection"))

df_raw <- df_raw %>% mutate(Fibrosis_Score = case_when(V19 == "Ishak 0-4; No to moderate fibrosis; METAVIR F0-F3; Batt-Ludwig 0-3" ~ "Ishak 0-4; No to moderate fibrosis",
                                                       V19 == "Ishak 5-6; Advanced/severe fibrosis; METAVIR F4; Batt-Ludwig 4; Cirrhosis" ~ "Ishak 5-6; Advanced/severe fibrosis",
                                                       TRUE ~ "Unknown"))

df_raw$AFP_status = ifelse(df_raw$V20 == "Positive/elevated", "Positive", "Negative")

df_raw$primary_noninvasive = ifelse(df_raw$V24 == "100", "Yes", "No")

df_raw$tumor_size = as.numeric(df_raw$V27) 

df_raw$OS_STATUS = ifelse(df_raw$V32 == "Dead",1,0)
df_raw$OS_MONTHS = as.numeric(df_raw$V33)


df_raw$first_malignant_primary = ifelse(df_raw$V34 == "Yes", "Yes", "No")

df_raw$num_malignant = as.numeric(df_raw$V35)

df_raw$Age = apply(df_raw,1,function(x){
  age_str = x["V37"]
  age_str2 = str_split(age_str," ")[[1]][1]
  if (age_str2 == "90+"){
    return(90)
  }else{
    age = as.numeric(str_split(age_str," ")[[1]][1])
    return(age)
  }
})

df_raw$Married = ifelse(df_raw$V38 == "Married (including common law)", "Yes", "No")

df_raw$adjusted_household_income = ifelse(df_raw$V39 %in% c("< $40,000", "$40,000 - $44,999", "$45,000 - $49,999", 
                                              "$50,000 - $54,999", "$55,000 - $59,999", "$60,000 - $64,999",
                                              "$65,000 - $69,999", "$70,000 - $74,999", "$75,000 - $79,999"), "<$80,000", "$80,000+")

df_raw <- df_raw[,22:42]

saveRDS(df_raw,file = "data/SEER_HCC_2018_2021_raw.rds")








