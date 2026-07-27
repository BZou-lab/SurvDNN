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
library(BART)

##### 0. Args Setting #####
use_condaenv("/nas/longleaf/home/{onyen}/.local/share/r-miniconda/envs/deepsurv/bin/python")
args1 = commandArgs(TRUE)
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) *2399
Fold_Name = "/work/users/{o}/{n}/{onyen}/Permfit_Sim/Realdata4/"
dir.create(paste(Fold_Name,Folder_args,sep = ""))
##### 1. Load RealData Set #####
Paindata_Use = readRDS("/work/users/{o}/{n}/{onyen}/Permfit_Sim/Realdata4/realdata4_METABRIC.rds")


#Paindata_Use$GSTATUS = as.numeric(Paindata_Use$GSTATUS)
#Paindata_Use$PTIME = Paindata_Use$PTIME - Paindata_Use$LOS

#Paindata_Use = Paindata_Use %>% filter(PTIME > 0)

PainData = importDnnetSurv(x = Paindata_Use[,3:71],
                           y = Paindata_Use$OS_MONTHS,
                           e = Paindata_Use$OS_STATUS)
use_seed = 256 +start.seed
set.seed(use_seed)

n_train = floor(nrow(Paindata_Use)*0.8)
train_id = sample(nrow(Paindata_Use))[1:n_train]
valid_id = which(c(1:nrow(Paindata_Use)) %in% train_id == F)

trainN = paste(Fold_Name,Folder_args,"/traindata","_",use_seed,".rds",sep = "")
validN = paste(Fold_Name,Folder_args,"/validdata","_",use_seed,".rds",sep = "")

if (trainN %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  PainData_train = importDnnetSurv(x = Paindata_Use[train_id,3:71],
                                   y = Paindata_Use[train_id,"OS_MONTHS"],
                                   e = Paindata_Use[train_id,"OS_STATUS"])
  
  PainData_valid = importDnnetSurv(x = Paindata_Use[valid_id,3:71],
                                   y = Paindata_Use[valid_id,"OS_MONTHS"],
                                   e = Paindata_Use[valid_id,"OS_STATUS"])
  
  saveRDS(Paindata_Use[train_id,],file = paste(Fold_Name,Folder_args,"/traindata","_",use_seed,".rds",sep = ""))
  saveRDS(Paindata_Use[valid_id,],file = paste(Fold_Name,Folder_args,"/validdata","_",use_seed,".rds",sep = ""))
}else{
  Ptrain = readRDS(trainN)
  Pvalid = readRDS(validN)
  PainData_train = importDnnetSurv(x = Ptrain[,3:71],
                                   y = Ptrain[,"OS_MONTHS"],
                                   e = Ptrain[,"OS_STATUS"])
  PainData_valid = importDnnetSurv(x = Pvalid[,3:71],
                                   y = Pvalid[,"OS_MONTHS"],
                                   e = Pvalid[,"OS_STATUS"])
}

######### Load Function ############

source("~/R-Code/source_permfit_survival_20230610.R")

# Perm DNNet 1 
FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet1","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dnnet1 = permfit_survival(train = PainData_train,n_perm =100,method = "dnnet",
                                 k_fold = 5,
                                 pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                     c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM")),
                                 n.hidden = c(50, 40, 30, 20),family = "coxph",
                                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                                 n.epoch = 1000, learning.rate.adaptive = "adam",
                                 plot = FALSE) %>% try()
  saveRDS(Perm_dnnet1,file = FileN)
  rm(Perm_dnnet1)
}
print(FileN)
print(Sys.time())
# Perm DNNet
FileN = paste(Fold_Name,Folder_args,"/PermFit_dnnet","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                 n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)
  Perm_dnnet = permfit_survival(train = PainData_train,n_perm =100,method = "ensemble_dnnet",
                                k_fold = 5,
                                pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                    c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM")),
                                n.ensemble = 100, esCtrl = esCtrl) %>% try()
  saveRDS(Perm_dnnet,file = FileN)
  rm(Perm_dnnet)
}
print(FileN)
print(Sys.time())
# Perm Xgboost
FileN = paste(Fold_Name,Folder_args,"/PermFit_xgb","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_xgb = permfit_survival(train = PainData_train,n_perm =100,method = "Xgboost",
                              k_fold = 5,
                              pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                  c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM")),
                              nrounds = 50) %>% try()
  saveRDS(Perm_xgb,file = FileN)
  rm(Perm_xgb)
}
print(FileN)
print(Sys.time())
# Perm RSF
FileN = paste(Fold_Name,Folder_args,"/PermFit_rsf","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_rsf = permfit_survival(train = PainData_train,n_perm =100,method = "random_forest",
                              k_fold = 5,
                              pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                  c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM")),
                              ntrees = 500) %>% try()
  saveRDS(Perm_rsf,file = FileN)
  rm(Perm_rsf)
}
print(FileN)
print(Sys.time())
# Perm COX
FileN = paste(Fold_Name,Folder_args,"/PermFit_cox","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_cox = permfit_survival(train = PainData_train,n_perm =100,method = "survival_cox",
                              k_fold = 5,
                              pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                  c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM"))) %>% try()
  saveRDS(Perm_cox,file = FileN)
  rm(Perm_cox)
}
print(FileN)
print(Sys.time())
# Perm AFT
FileN = paste(Fold_Name,Folder_args,"/PermFit_aft","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_aft = permfit_survival(train = PainData_train,n_perm =100,method = "survival_aft",
                              k_fold = 5,
                              pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                  c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM"))) %>% try()
  saveRDS(Perm_aft,file = FileN)
  rm(Perm_aft)
}
print(FileN)
print(Sys.time())
# Perm deephit
FileN = paste(Fold_Name,Folder_args,"/PermFit_dhit","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dhit = permfit_survival(train = PainData_train,n_perm =100,method = "DeepHit",
                               k_fold = 5,
                               pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                   c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM")),
                               activation = "relu",
                               frac = 0.2,early_stopping = T,
                               num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                               batch_size = 50) %>% try()
  saveRDS(Perm_dhit,file = FileN)
  rm(Perm_dhit)
}
print(FileN)
print(Sys.time())
# Perm deepsurv
FileN = paste(Fold_Name,Folder_args,"/PermFit_dsurv","_",use_seed,".rds",sep = "")

if (FileN%in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
  Perm_dsurv = permfit_survival(train = PainData_train,n_perm =100,method = "DeepSurv",
                                k_fold = 5,
                                pathway_list = list(c("CELLULARITY_M","CELLULARITY_H"),
                                                    c("CLAUDIN_SUBTYPE_AGGRESIVE","CLAUDIN_SUBTYPE_LUM")),
                                activation = "relu",
                                frac = 0.2,early_stopping = T,
                                num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                batch_size = 50) %>% try()
  saveRDS(Perm_dsurv,file = FileN)
  rm(Perm_dsurv)
}
print(FileN)
print(Sys.time())
# Perm SVM
#FileN = paste(Fold_Name,Folder_args,"/PermFit_svm","_",use_seed,".rds",sep = "")

#if (FileN%in% 
#list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F) {
#Perm_svm = permfit_survival(train = PainData_train,n_perm =100,method = "Survival_SVM",
#k_fold = 5,
# pathway_list = list(c("ABO_A","ABO_B","ABO_AB"),
#                     c("ABO_DON_A","ABO_DON_B","ABO_DON_AB"),
#                     c("COD_CAD_DON_A","COD_CAD_DON_C","COD_CAD_DON_O"),
#                     c("DIAG_IPF","DIAG_COPD")),
#  opt.meth = "quadprog",kernel = "add_kernel",
#  gamma.mu = 0.1,type = "vanbelle2",diff.meth = "makediff3") %>% try()
#saveRDS(Perm_svm,file = FileN)
#rm(Perm_svm)
#}
#print(FileN)
#print(Sys.time())