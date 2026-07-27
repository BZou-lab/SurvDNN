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

##### 0. Args Setting #####
use_condaenv("/nas/longleaf/home/{onyen}/.local/share/r-miniconda/envs/deepsurv/bin/python")

args1 = commandArgs(TRUE)
Folder_args = as.numeric(args1[1])
start.seed = as.numeric(args1[1]) * 1999

# Working Directory
Fold_Name = "/work/users/{o}/{n}/{onyen}/Permfit_Sim/formal_newpermfits3_lowcor_1000_5/"
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

DataN_Scov = paste(Fold_Name,Folder_args,"/Data_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")


Full_Scov = readRDS(DataN_Scov)


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


##### 3. Load Functions #####

##### 4. Permfit Survival #####

source("/nas/longleaf/home/{onyen}/R-Code/source_permfit_survival_20230610.R")

p_thershold = 0.05

##### DNNet1 #####

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dnnet1_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "dnnet",model.type = "survival",
                            object = Train_Scov,
                            n.hidden = c(50, 40, 30, 20),family = "coxph",
                            l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                            n.epoch = 1000, learning.rate.adaptive = "adam", 
                            plot = FALSE) %>% try()
  
  selected_model_S = mod_permfit(method = "dnnet",model.type = "survival",
                                 object = Selected_Train_Scov,
                                 n.hidden = c(50, 40, 30, 20),family = "coxph",
                                 l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                                 n.epoch = 1000, learning.rate.adaptive = "adam", 
                                 plot = FALSE) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "dnnet",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "dnnet",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Cdnnet1_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Cdnnet1_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Cdnnet1_org_S = NA
  Cdnnet1_selected_S = NA
}


print(Sys.time())

##### DNNet #####

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dnnet_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
               l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
               n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                            object = Train_Scov,
                            n.ensemble = 100, esCtrl = esCtrl) %>% try()
  
  selected_model_S = mod_permfit(method = "ensemble_dnnet",model.type = "survival",
                                 object = Selected_Train_Scov,
                                 n.ensemble = 100, esCtrl = esCtrl) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "ensemble_dnnet",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "ensemble_dnnet",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Cdnnet_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Cdnnet_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Cdnnet_org_S = NA
  Cdnnet_selected_S = NA
}


print(Sys.time())

##### XGBOOST #####

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_xgb_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "Xgboost",model.type = "survival",
                            object = Train_Scov,
                            nrounds = 50) %>% try()
  
  selected_model_S = mod_permfit(method = "Xgboost",model.type = "survival",
                                 object = Selected_Train_Scov,
                                 nrounds = 50) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "Xgboost",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "Xgboost",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Cxgb_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Cxgb_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Cxgb_org_S = NA
  Cxgb_selected_S = NA
}

print(Sys.time())

##### RSF #####

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_rsf_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "random_forest",model.type = "survival",
                            object = Train_Scov,
                            ntrees = 500) %>% try()
  
  selected_model_S = mod_permfit(method = "random_forest",model.type = "survival",
                                 object = Selected_Train_Scov,
                                 ntrees = 500) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "random_forest",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "random_forest",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Crsf_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Crsf_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Crsf_org_S = NA
  Crsf_selected_S = NA
}


print(Sys.time())

##### COX #####

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_cox_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "survival_cox",model.type = "survival",
                            object = Train_Scov) %>% try()
  
  selected_model_S = mod_permfit(method = "survival_cox",model.type = "survival",
                                 object = Selected_Train_Scov) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "survival_cox",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "survival_cox",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Ccox_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Ccox_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Ccox_org_S = NA
  Ccox_selected_S = NA
}


print(Sys.time())

##### AFT #####

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_aft_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "survival_aft",model.type = "survival",
                            object = Train_Scov) %>% try()
  
  selected_model_S = mod_permfit(method = "survival_aft",model.type = "survival",
                                 object = Selected_Train_Scov) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "survival_aft",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "survival_aft",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Caft_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Caft_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Caft_org_S = NA
  Caft_selected_S = NA
}


print(Sys.time())

##### DeepSurv #####


FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dsurv_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "DeepSurv",model.type = "survival",
                            object = Train_Scov,activation = "relu",
                            frac = 0.2,early_stopping = T,
                            num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                            batch_size = 50) %>% try()
  
  selected_model_S = mod_permfit(method = "DeepSurv",model.type = "survival",
                                 object = Selected_Train_Scov,activation = "relu",
                                 frac = 0.2,early_stopping = T,
                                 num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                 batch_size = 50) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "DeepSurv",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "DeepSurv",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Cdsurv_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Cdsurv_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Cdsurv_org_S = NA
  Cdsurv_selected_S = NA
}


print(Sys.time())

##### DeepHIT #####


FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dhit_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

result_S = readRDS(file = FileN_Scov) %>% try()

if(class(result_S) != "try-error"){
  imp_S = result_S@importance
  imp_block_S = result_S@block_importance
  sig_var_S = imp_S$var_name[which(imp_S$importance_pval_x < p_thershold)]
  
  if(length(sig_var_S) < 2){
    arrange_S = imp_S %>% arrange(importance_pval_x)
    sig_var_S = arrange_S$var_name[1:6]
  }
  
  Selected_Train_Scov = importDnnetSurv(x = Full_Scov[1:N_train,sig_var_S],
                                        y = Full_Scov$eventtime1[1:N_train],
                                        e = Full_Scov$status1[1:N_train])
  Selected_Valid_Scov = importDnnetSurv(x = Full_Scov[(N_train+1):N,sig_var_S],
                                        y = Full_Scov$eventtime1[(N_train+1):N],
                                        e = Full_Scov$status1[(N_train+1):N])
  org_model_S = mod_permfit(method = "DeepHit",model.type = "survival",
                            object = Train_Scov,activation = "relu",
                            frac = 0.2,early_stopping = T,
                            num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                            batch_size = 50) %>% try()
  
  selected_model_S = mod_permfit(method = "DeepHit",model.type = "survival",
                                 object = Selected_Train_Scov,activation = "relu",
                                 frac = 0.2,early_stopping = T,
                                 num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                 batch_size = 50) %>% try()
  
  pred_org_S = predict_mod_permfit(mod = org_model_S,
                                   method = "DeepHit",
                                   model.type = "survival",
                                   object = Valid_Scov) %>% try()
  
  pred_selected_S = predict_mod_permfit(mod = selected_model_S,
                                        method = "DeepHit",
                                        model.type = "survival",
                                        object = Selected_Valid_Scov) %>% try()
  
  Cdhit_org_S = try(1 - rcorr.cens(pred_org_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
  Cdhit_selected_S = try(1 - rcorr.cens(pred_selected_S[[1]],Surv(Valid_Scov@y,Valid_Scov@e))[1])
}else {
  Cdhit_org_S = NA
  Cdhit_selected_S = NA
}



print(Sys.time())

c_all_S = list(Ccox_org_S,Ccox_selected_S,
               Caft_org_S,Caft_selected_S,
               Cdnnet1_org_S,Cdnnet1_selected_S,
               Cdnnet_org_S,Cdnnet_selected_S,
               Crsf_org_S,Crsf_selected_S,
               Cxgb_org_S,Cxgb_selected_S,
               #Csvm_org,Csvm_selected,
               Cdsurv_org_S,Cdsurv_selected_S,
               Cdhit_org_S,Cdhit_selected_S)


saveRDS(c_all_S,file = paste(Fold_Name,Folder_args,"/Summary_prediction_005_S.rds",sep = ""))










