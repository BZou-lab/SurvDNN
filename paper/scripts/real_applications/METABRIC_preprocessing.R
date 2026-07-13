# Script to clean and preprocess METABRIC data for survival analysis

rm(list = ls())

library(data.table)
library(tidyverse)
library(dplyr)
library(survival)
library(randomForestSRC)
library(glmnet)

meta_clin_p = as.data.frame(fread("data/data_clinical_patient.txt"))
colnames(meta_clin_p) = meta_clin_p[4,]
meta_clin_p = meta_clin_p[-c(1:4),]
meta_clin_s = as.data.frame(fread("data/data_clinical_sample.txt"))
colnames(meta_clin_s) = meta_clin_s[4,]
meta_clin_s = meta_clin_s[-c(1:4),]

meta_mut = as.data.frame(fread("data/data_mutations.txt"))
mut_DF = table(meta_mut$Hugo_Symbol) %>% as.data.frame() %>% filter(Freq >= 100)
High_Mut_Genes = mut_DF$Var1 %>% as.character()

meta_clin_ps = meta_clin_p %>% 
  filter(LYMPH_NODES_EXAMINED_POSITIVE != "") %>%
  filter(NPI != "") %>%
  filter(CELLULARITY != "") %>%
  filter(CHEMOTHERAPY != "") %>%
  filter(ER_IHC != "") %>%
  filter(HER2_SNP6 != "" & HER2_SNP6 != "UNDEF") %>%
  filter(HORMONE_THERAPY != "") %>%
  filter(INFERRED_MENOPAUSAL_STATE != "") %>%
  filter(SEX != "") %>%
  filter(AGE_AT_DIAGNOSIS != "") %>%
  filter(OS_STATUS != "") %>%
  filter(OS_MONTHS != "") %>%
  filter(LATERALITY != "") %>%
  filter(RADIO_THERAPY != "") %>%
  filter(HISTOLOGICAL_SUBTYPE != "") %>%
  filter(BREAST_SURGERY != "") %>%
  filter(RFS_STATUS != "") %>%
  select(c(PATIENT_ID,OS_STATUS,OS_MONTHS,LYMPH_NODES_EXAMINED_POSITIVE,NPI,CELLULARITY,
           CHEMOTHERAPY,ER_IHC,HER2_SNP6,HORMONE_THERAPY,INFERRED_MENOPAUSAL_STATE,AGE_AT_DIAGNOSIS,CLAUDIN_SUBTYPE,LATERALITY,RADIO_THERAPY,
           HISTOLOGICAL_SUBTYPE,BREAST_SURGERY,RFS_STATUS))

meta_clin_ps$OS_STATUS = ifelse(meta_clin_ps$OS_STATUS == "1:DECEASED",1,0)
meta_clin_ps$OS_MONTHS = as.numeric(meta_clin_ps$OS_MONTHS)
hist(as.numeric(meta_clin_ps$LYMPH_NODES_EXAMINED_POSITIVE))
meta_clin_ps$LYMPH_NODES_EXAMINED_POSITIVE = as.numeric(meta_clin_ps$LYMPH_NODES_EXAMINED_POSITIVE)
hist(as.numeric(meta_clin_ps$NPI))
meta_clin_ps$NPI = as.numeric(meta_clin_ps$NPI)
table(meta_clin_ps$CELLULARITY,useNA = "always")
meta_clin_ps$CELLULARITY_M = ifelse(meta_clin_ps$CELLULARITY == "Moderate",1,0)
meta_clin_ps$CELLULARITY_H = ifelse(meta_clin_ps$CELLULARITY == "High",1,0)
table(meta_clin_ps$CHEMOTHERAPY)
meta_clin_ps$CHEMOTHERAPY_Y = ifelse(meta_clin_ps$CHEMOTHERAPY == "YES",1,0)
table(meta_clin_ps$ER_IHC)
meta_clin_ps$ER_IHC_P = ifelse(meta_clin_ps$ER_IHC == "Positve",1,0)
table(meta_clin_ps$HER2_SNP6)
meta_clin_ps$HER2_SNP6_ABNORMAL = ifelse(meta_clin_ps$HER2_SNP6 !="NEUTRAL",1,0)
table(meta_clin_ps$HORMONE_THERAPY)
meta_clin_ps$HORMONE_THERAPY_Y = ifelse(meta_clin_ps$HORMONE_THERAPY == "YES",1,0)
table(meta_clin_ps$INFERRED_MENOPAUSAL_STATE)
meta_clin_ps$INFERRED_MENOPAUSAL_STATE_PRE = ifelse(meta_clin_ps$INFERRED_MENOPAUSAL_STATE == "Pre",1,0)
meta_clin_ps$AGE_AT_DIAGNOSIS = as.numeric(meta_clin_ps$AGE_AT_DIAGNOSIS)
table(meta_clin_ps$CLAUDIN_SUBTYPE)
meta_clin_ps$CLAUDIN_SUBTYPE_AGGRESIVE = ifelse(meta_clin_ps$CLAUDIN_SUBTYPE == "Basal" | 
                                                  meta_clin_ps$CLAUDIN_SUBTYPE == "claudin-low" | 
                                                  meta_clin_ps$CLAUDIN_SUBTYPE == "Her2" | 
                                                  meta_clin_ps$CLAUDIN_SUBTYPE == "NC",1,0)
meta_clin_ps$CLAUDIN_SUBTYPE_LUM = ifelse(meta_clin_ps$CLAUDIN_SUBTYPE == "LumA" | 
                                            meta_clin_ps$CLAUDIN_SUBTYPE == "LumB",1,0)
table(meta_clin_ps$LATERALITY) 
meta_clin_ps$LATERALITY_R = ifelse(meta_clin_ps$LATERALITY == "Right",1,0)
table(meta_clin_ps$RADIO_THERAPY)
meta_clin_ps$RADIO_THERAPY_Y = ifelse(meta_clin_ps$RADIO_THERAPY == "YES",1,0)
table(meta_clin_ps$HISTOLOGICAL_SUBTYPE)
meta_clin_ps$HISTOLOGICAL_SUBTYPE_COMMON = ifelse(meta_clin_ps$HISTOLOGICAL_SUBTYPE == "Ductal/NST" | 
                                                    meta_clin_ps$HISTOLOGICAL_SUBTYPE == "Lobular" |
                                                    meta_clin_ps$HISTOLOGICAL_SUBTYPE == "Tubular/ cribriform",1,0)
table(meta_clin_ps$BREAST_SURGERY)
meta_clin_ps$BREAST_SURGERY_M = ifelse(meta_clin_ps$BREAST_SURGERY == "MASTECTOMY",1,0)
table(meta_clin_ps$RFS_STATUS)
meta_clin_ps$RFS_STATUS_R = ifelse(meta_clin_ps$RFS_STATUS == "1:Recurred",1,0)

meta_clin_ps = meta_clin_ps %>%
  select(c(PATIENT_ID,OS_STATUS,OS_MONTHS,LYMPH_NODES_EXAMINED_POSITIVE,NPI,AGE_AT_DIAGNOSIS,
           CELLULARITY_M,CELLULARITY_H,CHEMOTHERAPY_Y,ER_IHC_P,HER2_SNP6_ABNORMAL,
           HORMONE_THERAPY_Y,INFERRED_MENOPAUSAL_STATE_PRE,CLAUDIN_SUBTYPE_AGGRESIVE,
           CLAUDIN_SUBTYPE_LUM,LATERALITY_R,RADIO_THERAPY_Y,HISTOLOGICAL_SUBTYPE_COMMON,
           BREAST_SURGERY_M,RFS_STATUS_R))

meta_clin_ss = meta_clin_s %>%
  filter(ER_STATUS != "") %>% 
  filter(!is.na(ER_STATUS)) %>%
  filter(HER2_STATUS != "") %>% 
  filter(PR_STATUS != "") %>% 
  filter(!is.na(TUMOR_SIZE)) %>% 
  #filter(!is.na(TUMOR_STAGE)) %>% 
  filter(!is.na(TMB_NONSYNONYMOUS)) %>%
  select(c(PATIENT_ID,ER_STATUS,HER2_STATUS,PR_STATUS,TUMOR_SIZE,TMB_NONSYNONYMOUS))

table(meta_clin_ss$ER_STATUS)
meta_clin_ss$ER_STATUS_P = ifelse(meta_clin_ss$ER_STATUS == "Positive",1,0)
table(meta_clin_ss$HER2_STATUS)
meta_clin_ss$HER2_STATUS_P = ifelse(meta_clin_ss$HER2_STATUS == "Positive",1,0)
table(meta_clin_ss$PR_STATUS)
meta_clin_ss$PR_STATUS_P = ifelse(meta_clin_ss$PR_STATUS == "Positive",1,0)
table(meta_clin_ss$TMB_NONSYNONYMOUS)
meta_clin_ss$TUMOR_SIZE = as.numeric(meta_clin_ss$TUMOR_SIZE)
meta_clin_ss$TMB_NONSYNONYMOUS = as.numeric(meta_clin_ss$TMB_NONSYNONYMOUS)

meta_clin_ss = meta_clin_ss %>%
  select(c(PATIENT_ID,ER_STATUS_P,HER2_STATUS_P,PR_STATUS_P,TUMOR_SIZE,TMB_NONSYNONYMOUS))

meta_survival = inner_join(meta_clin_ps,meta_clin_ss,by = "PATIENT_ID")

meta_survival_mut = cbind(meta_survival,as.data.frame(matrix(NA,nrow = nrow(meta_survival),
                                                             ncol = length(High_Mut_Genes))))

colnames(meta_survival_mut)[(ncol(meta_survival_mut)-length(High_Mut_Genes)+1):ncol(meta_survival_mut)] = High_Mut_Genes

for (i in 1:nrow(meta_survival_mut)){
  sampleID = meta_survival_mut[i,"PATIENT_ID"]
  sample_mut = meta_mut %>% filter(Tumor_Sample_Barcode == sampleID)
  if (nrow(sample_mut) == 0){
    meta_survival_mut[i,(ncol(meta_survival_mut)-length(High_Mut_Genes)+1):ncol(meta_survival_mut)] = 0
  }else{
    sample_mut_gene = unique(sample_mut$Hugo_Symbol)
    meta_survival_mut[i,High_Mut_Genes[which(High_Mut_Genes %in% sample_mut_gene)]] = 1
    meta_survival_mut[i,High_Mut_Genes[which(!High_Mut_Genes %in% sample_mut_gene)]] = 0
  }
}

meta_survival_mut = meta_survival_mut[,-1]

saveRDS(meta_survival_mut,file = "data/realdata4_METABRIC.rds")