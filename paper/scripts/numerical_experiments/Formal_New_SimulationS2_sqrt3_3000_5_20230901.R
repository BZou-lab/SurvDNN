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
Fold_Name = "/work/users/{o}/{n}/{onyen}/Permfit_Sim/formal_newpermfits2_lowcor_3000_5_sqrt3/"
dir.create(paste(Fold_Name,Folder_args,sep = ""))

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

##### 1.  Simulate Covariate with covarience ######

# Data's Name

DataN_Scov = paste(Fold_Name,Folder_args,"/Data_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (DataN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == T){
  Full_Scov = readRDS(DataN_Scov)
}else{
  rho = 0.5
  
  X_Scov = data.frame(mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho),
                      mvrnorm(n=N,mu=rep(0,p_block),
                              Sigma = (1-rho)*diag(p_block)+rho))
  colnames(X_Scov) = paste("X",1:(4*p_block),sep = "")
  
  for (z in 1:p_block) {
    if (z==1){
      Z_Scov = data.frame(rbinom(N,1,0.4))
    }else{
      Z_random = data.frame(rbinom(N,1,0.4))
      Z_Scov = cbind(Z_Scov,Z_random)
    }
  }
  
  colnames(Z_Scov) = paste("Z",1:(p_block),sep = "")
  
  C_Scov = cbind(X_Scov,Z_Scov)
  
  # Significant variable: Z1,X1,X11,X21,X31
  
  sims_beta = c(0.2,0.4,0.6,0.8,1,rep(0,p-5))/sqrt(3)
  
  names(sims_beta) = c("X1",
                       paste("X",p_block+1,sep = ""),
                       paste("X",2*p_block+1,sep = ""),
                       paste("X",3*p_block+1,sep = ""),
                       "Z1",paste("X",(1:(4*p_block))[-c(1,p_block+1,2*p_block+1,3*p_block+1)],sep = ""),
                       paste("Z",2:p_block,sep = ""))
  
  # Generate Non-Linear Term
  
  #C_Scov[,"X1Z2"] = C_Scov[,"X1"]*C_Scov[,"Z2"]
  #C_Scov[,paste("X",p_block+1,"Square",sep = "")] = (C_Scov[,paste("X",p_block+1,sep = "")])^2
  #C_Scov[,paste("X",2*p_block+1,"X",3*p_block+1,sep = "")] = C_Scov[,paste("X",2*p_block+1,sep = "")]*C_Scov[,paste("X",3*p_block+1,sep = "")]
  
  # Generate Survival Time
  
  Y_Scov= simsurv(dist = "gompertz",lambdas = lambda.t,gammas = alpha.t,
                  betas = sims_beta,
                  x= C_Scov)
  
  C_Scov$id=1:nrow(C_Scov)
  Full_Scov = left_join(Y_Scov,C_Scov,by = "id")
  
  # Generate Censor Distribution
  
  # Censor Distribution
  
  theta = 3.5
  
  set.seed(use_seed)
  
  Full_Scov$censor_time = rweibull(N,shape = alpha.c,scale = theta)
  
  Full_Scov$eventtime1 = apply(Full_Scov,1,function(x){
    Cens = x["censor_time"]
    Event = x["eventtime"]
    return(min(c(Cens,Event)))
  })
  
  Full_Scov$status1 = apply(Full_Scov,1,function(x){
    Cens = x["censor_time"]
    Event = x["eventtime"]
    if (Event >= Cens){
      return(0)
    } else {
      return(1)
    }
  })
  
  saveRDS(Full_Scov,file = DataN_Scov)
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


##### 3. Load Functions #####

##### 4. Permfit Survival #####

source("/nas/longleaf/home/{onyen}/R-Code/source_permfit_survival_20230610.R")

##### 5. Run Simulation, simulation result

# PermFit - DNNet 1

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dnnet1_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_dnnet1_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "dnnet",
                                      k_fold = 5,
                                      n.hidden = c(50, 40, 30, 20),family = "coxph",
                                      l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
                                      n.epoch = 1000, learning.rate.adaptive = "adam", 
                                      plot = FALSE) %>% try()
  saveRDS(Perm_dnnet1_Scov,file = FileN_Scov)
  rm(Perm_dnnet1_Scov)
}



print(Sys.time())

# PermFit - DNNet

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dnnet_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")
esCtrl <- list(n.hidden = c(50, 40, 30, 20), activate = "relu",
               l1.reg = 10**-4, early.stop.det = 1000, n.batch = 50,
               n.epoch = 1000, learning.rate.adaptive = "adam", plot = FALSE)

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_dnnet_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "ensemble_dnnet",
                                     k_fold = 5,
                                     n.ensemble = 100, esCtrl = esCtrl) %>% try()
  saveRDS(Perm_dnnet_Scov,file = FileN_Scov)
  rm(Perm_dnnet_Scov)
}



print(Sys.time())

# PermFit - Xgboost
FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_xgb_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_xgb_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "Xgboost",
                                   k_fold = 5,nrounds = 50) %>% try()
  saveRDS(Perm_xgb_Scov,file = FileN_Scov)
  rm(Perm_xgb_Scov)
}


print(Sys.time())

# PermFit - Random Survival Forest

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_rsf_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_rsf_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "random_forest",
                                   k_fold = 5,ntrees = 500) %>% try()
  saveRDS(Perm_rsf_Scov,file = FileN_Scov)
  rm(Perm_rsf_Scov)
}


print(Sys.time())

# PermFit - SVM

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_svm_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_svm_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "Survival_SVM",
                                   k_fold = 5,opt.meth = "quadprog",kernel = "lin_kernel",
                                   gamma.mu = 0.1,type = "vanbelle2",diff.meth = "makediff3") %>% try()
  saveRDS(Perm_svm_Scov,file = FileN_Scov)
  rm(Perm_svm_Scov)
}

print(Sys.time())

# PermFit - Cox

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_cox_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")


if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_cox_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "survival_cox",
                                   k_fold = 5) %>% try()
  saveRDS(Perm_cox_Scov,file = FileN_Scov)
  rm(Perm_cox_Scov)
}

print(Sys.time())

# PermFit - AFT 

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_aft_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_aft_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "survival_aft",
                                   k_fold = 5) %>% try()
  saveRDS(Perm_aft_Scov,file = FileN_Scov)
  rm(Perm_aft_Scov)
}


# PermFit - DeepHit

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dhit_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_deephit_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "DeepHit",
                                       k_fold = 5,activation = "relu",
                                       frac = 0.2,early_stopping = T,
                                       num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                       batch_size = 50) %>% try()
  saveRDS(Perm_deephit_Scov,file = FileN_Scov)
  rm(Perm_deephit_Scov)
}

print(Sys.time())

# PermFit - DeepSurv 

FileN_Scov = paste(Fold_Name,Folder_args,"/PermFit_dsurv_Scov","_",
                   sim_num,"_",use_seed,"_",p_block,".rds",sep = "")

if (FileN_Scov %in% 
    list.files(paste(Fold_Name,Folder_args,sep = ""),full.names = T) == F){
  Perm_deepsurv_Scov = permfit_survival(train = Train_Scov,n_perm =100,method = "DeepSurv",
                                        k_fold = 5,activation = "relu",
                                        frac = 0.2,early_stopping = T,
                                        num_nodes = c(50L, 40L, 30L,20L),epochs = 1000,
                                        batch_size = 50) %>% try()
  saveRDS(Perm_deepsurv_Scov,file = FileN_Scov)
  rm(Perm_deepsurv_Scov)
}

print(Sys.time())


