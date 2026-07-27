rm(list = ls())
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
library(shapr)
library(SHAPforxgboost)
library(praznik)

##### 0. Args Setting #####
use_condaenv("/nas/longleaf/home/{onyen}/.local/share/r-miniconda/envs/deepsurv/bin/python")

args1 = commandArgs(TRUE)
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) * 1999

# Working Directory
Fold_Name = "/work/users/{o}/{n}/{onyen}/Permfit_Sim/formal_newpermfits2_lowcor_3000_sqrt3/"
# dir.create(paste(Fold_Name,Folder_args,sep = ""))

# N = 1200
# p_block = 20

N = args1[2] %>% as.numeric()
p_block = args1[3] %>% as.numeric()

N_train = N*5/6

p = 5*p_block

# Survival Simulation Parameter
alpha.t = 3
lambda.t = 0.00002
alpha.c = 1.5

# Simulation Seed

sim_num = 1
use_seed = sim_num+start.seed
set.seed(use_seed)

##### 2. Load Functions #####

source("/nas/longleaf/home/{onyen}/R-Code/source_permfit_survival_20230610.R")

##### 3. Load Data ######

DataN_Ncov = paste(Fold_Name,Folder_args,"/Data_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (DataN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == T){
  Full_Ncov = readRDS(DataN_Ncov)
}

DataN_Scov = paste(Fold_Name,Folder_args,"/Data_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (DataN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == T){
  Full_Scov = readRDS(DataN_Scov)
}

Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Scov$eventtime1[1:N_train],
                             e = Full_Scov$status1[1:N_train])
Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Scov$eventtime1[(N_train+1):N],
                             e = Full_Scov$status1[(N_train+1):N])

Train_Ncov = importDnnetSurv(x = Full_Ncov[1:N_train,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Ncov$eventtime1[1:N_train],
                             e = Full_Ncov$status1[1:N_train])
Valid_Ncov = importDnnetSurv(x = Full_Ncov[(N_train+1):N,
                                           c(paste("X",1:(4*p_block),sep = ""),
                                             paste("Z",1:(p_block),sep = ""))],
                             y = Full_Ncov$eventtime1[(N_train+1):N],
                             e = Full_Ncov$status1[(N_train+1):N])

#### 4. Original Variable Selection Method #####

# Original - XGBoost

FileN_Scov = paste(Fold_Name,Folder_args,"/Org_xgb_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Scov2 = paste(Fold_Name,Folder_args,"/Org_xgbshap_Scov","_",
                    sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov = paste(Fold_Name,Folder_args,"/Org_xgb_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov2 = paste(Fold_Name,Folder_args,"/Org_xgbshap_Ncov","_",
                    sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_xgb_Scov = mod_permfit(method = "Xgboost",model.type = "survival",
                             object = Train_Scov,nrounds = 50) %>% try()
  # XGboost Importance
  xgb_imp_Scov = xgb.importance(model = Org_xgb_Scov)
  saveRDS(xgb_imp_Scov,file = FileN_Scov)
  # SHAP-XGboost
  # To return the SHAP values and ranked features by mean|SHAP|
  shap_values_Scov <- shap.values(xgb_model = Org_xgb_Scov, X_train = Train_Scov@x) %>% try()
  # The ranked features by mean |SHAP|
  # shap_values$mean_shap_score
  saveRDS(shap_values_Scov,file = FileN_Scov2)
  rm(Org_xgb_Scov)
  print(Sys.time())
}

if (FileN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_xgb_Ncov = mod_permfit(method = "Xgboost",model.type = "survival",
                             object = Train_Ncov,nrounds = 50) %>% try()
  # XGboost Importance
  xgb_imp_Ncov = xgb.importance(model = Org_xgb_Ncov)
  saveRDS(xgb_imp_Ncov,file = FileN_Ncov)
  # SHAP-XGboost
  # To return the SHAP values and ranked features by mean|SHAP|
  shap_values_Ncov <- shap.values(xgb_model = Org_xgb_Ncov, X_train = Train_Ncov@x) %>% try()
  # The ranked features by mean |SHAP|
  # shap_values$mean_shap_score
  saveRDS(shap_values_Ncov,file = FileN_Ncov2)
  rm(Org_xgb_Ncov)
  print(Sys.time())
}

print(Sys.time())

# Original - RSF 

FileN_Scov = paste(Fold_Name,Folder_args,"/Org_rsf_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov = paste(Fold_Name,Folder_args,"/Org_rsf_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_rsf_Scov = mod_permfit(method = "random_forest",model.type = "survival",
                             object = Train_Scov,ntrees = 500) %>% try()
  # RSF Importance
  subrsf = subsample(Org_rsf_Scov) %>% try()
  Sub_rsf_Scov = extract.subsample(subrsf) %>% try()
  saveRDS(Sub_rsf_Scov,file = FileN_Scov)
  rm(Org_rsf_Scov)
}

if (FileN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_rsf_Ncov = mod_permfit(method = "random_forest",model.type = "survival",
                             object = Train_Ncov,ntrees = 500) %>% try()
  # RSF Importance
  subrsf = subsample(Org_rsf_Ncov) %>% try()
  Sub_rsf_Ncov = extract.subsample(subrsf) %>% try()
  saveRDS(Sub_rsf_Ncov,file = FileN_Ncov)
  rm(Org_rsf_Ncov)
}
print(Sys.time())
# Original - COX 

FileN_Scov = paste(Fold_Name,Folder_args,"/Org_cox_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov = paste(Fold_Name,Folder_args,"/Org_cox_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_cox_Scov = mod_permfit(method = "survival_cox",model.type = "survival",
                             object = Train_Scov) %>% try()
  saveRDS(Org_cox_Scov,file = FileN_Scov)
  rm(Org_cox_Scov)
}

if (FileN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_cox_Ncov = mod_permfit(method = "survival_cox",model.type = "survival",
                             object = Train_Ncov) %>% try()
  saveRDS(Org_cox_Ncov,file = FileN_Ncov)
  rm(Org_cox_Ncov)
}
print(Sys.time())
# Original - LASSO-COX

FileN_Scov = paste(Fold_Name,Folder_args,"/Org_coxlasso_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov = paste(Fold_Name,Folder_args,"/Org_coxlasso_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_coxlasso_Scov = mod_permfit(method = "lasso",model.type = "survival",
                                  object = Train_Scov) %>% try()
  saveRDS(Org_coxlasso_Scov,file = FileN_Scov)
  rm(Org_coxlasso_Scov)
}

if (FileN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_coxlasso_Ncov = mod_permfit(method = "lasso",model.type = "survival",
                                  object = Train_Ncov) %>% try()
  saveRDS(Org_coxlasso_Ncov,file = FileN_Ncov)
  rm(Org_coxlasso_Ncov)
}
print(Sys.time())
# Original - AFT

FileN_Scov = paste(Fold_Name,Folder_args,"/Org_aft_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov = paste(Fold_Name,Folder_args,"/Org_aft_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_aft_Scov = mod_permfit(method = "survival_aft",model.type = "survival",
                             object = Train_Scov) %>% try()
  saveRDS(Org_aft_Scov,file = FileN_Scov)
  rm(Org_aft_Scov)
}

if (FileN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_aft_Ncov = mod_permfit(method = "survival_aft",model.type = "survival",
                             object = Train_Ncov) %>% try()
  saveRDS(Org_aft_Ncov,file = FileN_Ncov)
  rm(Org_aft_Ncov)
}
print(Sys.time())
# Original - Information Based

FileN_Scov = paste(Fold_Name,Folder_args,"/Org_JMIM_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
FileN_Ncov = paste(Fold_Name,Folder_args,"/Org_JMIM_Ncov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_cox_Scov = mod_permfit(method = "survival_cox",model.type = "survival",
                             object = Train_Scov) %>% try()
  Mar_Res_Scov = residuals(Org_cox_Scov,type = "martingale")
  q = max(min(N_train/3,10),2)
  Cut_MRes_Scov = cut(Mar_Res_Scov,breaks = 10)
  Train_Scov_Cat = as.data.frame(rep(NA,N_train))
  for (i in 1:(p-p_block)) {
    Cut_X = cut(Train_Scov@x[,i],breaks = 10) %>% as.data.frame()
    Train_Scov_Cat = cbind(Train_Scov_Cat,Cut_X)
  }
  for (i in (p-p_block+1):p) {
    Train_Scov_Cat = cbind(Train_Scov_Cat,as.factor(Train_Scov@x[,i]))
  }
  Train_Scov_Cat = Train_Scov_Cat[,-1]
  colnames(Train_Scov_Cat) = c(paste("X",1:(4*p_block),sep = ""),paste("Z",1:(p_block),sep = ""))
  JMIM_Scov = JMIM(Train_Scov_Cat,Cut_MRes_Scov,k=p)
  saveRDS(JMIM_Scov,file = FileN_Scov)
  rm(Org_cox_Scov)
  rm(JMIM_Scov)
}

if (FileN_Ncov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Org_cox_Ncov = mod_permfit(method = "survival_cox",model.type = "survival",
                             object = Train_Ncov) %>% try()
  Mar_Res_Ncov = residuals(Org_cox_Ncov,type = "martingale")
  q = max(min(N_train/3,10),2)
  Cut_MRes_Ncov = cut(Mar_Res_Ncov,breaks = 10)
  Train_Ncov_Cat = as.data.frame(rep(NA,N_train))
  for (i in 1:(p-p_block)) {
    Cut_X = cut(Train_Ncov@x[,i],breaks = 10) %>% as.data.frame()
    Train_Ncov_Cat = cbind(Train_Ncov_Cat,Cut_X)
  }
  for (i in (p-p_block+1):p) {
    Train_Ncov_Cat = cbind(Train_Ncov_Cat,as.factor(Train_Ncov@x[,i]))
  }
  Train_Ncov_Cat = Train_Ncov_Cat[,-1]
  colnames(Train_Ncov_Cat) = c(paste("X",1:(4*p_block),sep = ""),paste("Z",1:(p_block),sep = ""))
  JMIM_Ncov = JMIM(Train_Ncov_Cat,Cut_MRes_Ncov,k=p)
  saveRDS(JMIM_Ncov,file = FileN_Ncov)
  rm(Org_cox_Ncov)
  rm(JMIM_Ncov)
}
print(Sys.time())
