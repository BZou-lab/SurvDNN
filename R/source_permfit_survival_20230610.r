#if (!require("survival")){install.packages("survival")}else{library(survival)}
#library(survival)
#if (!require("tidyverse")){install.packages("tidyverse")}else{library(tidyverse)}
#library(tidyverse)
#if (!require("simsurv")){install.packages("simsurv")}else{library(simsurv)}
#library(simsurv)
#if (!require("MASS")){install.packages("MASS")}else{library(MASS)}
#library(MASS)
#if (!require("deepTL")){devtools::install_github("SkadiEye/deepTL")}else{library(deepTL)}
#library(deepTL)
#if (!require("survivalmodels")){install.packages("survivalmodels")}else{library(survivalmodels)}
#library(survivalmodels)
#if (!require("randomForestSRC")){install.packages("randomForestSRC")}else{library(randomForestSRC)}
#library(randomForestSRC)
#if (!require("survivalsvm")){install.packages("survivalsvm")}else{library(survivalsvm)}
#library(survivalsvm)
#if (!require("reticulate")){install.packages("reticulate")}else{library(reticulate)}
#library(reticulate)
#if (!require("Hmisc")){install.packages("Hmisc")}else{library(Hmisc)}
#library(Hmisc)
#if (!require("xgboost")){install.packages("xgboost")}else{library(xgboost)}
#library(xgboost)
#if (!require("pheatmap")){install.packages("pheatmap")}else{library(pheatmap)}
#library(pheatmap)

#' Title Check SurvDeepFIT's dependency, especially for DeepSurv and DeepHit
#'
#' @param package The name of the package to check dependencies for. Default is "SurvDeepFIT".
#' @param use_condaenv The name of the Conda environment to use, as specified by reticulate. Can be changed to a specific Conda environment.
#'
#' @return This function checks the dependencies for SurvDeepFIT, DeepSurv, DeepHit, and other related machine learning methods.
#' @export
#' @author Shiyu Wan
#'
#' @examples check_dependency()
check_dependency = function(package = "SurvDeepFIT",use_condaenv = "r-reticulate"){
  if (!require("survival")){install.packages("survival")}else{library(survival)}
  #library(survival)
  if (!require("tidyverse")){install.packages("tidyverse")}else{library(tidyverse)}
  #library(tidyverse)
  if (!require("simsurv")){install.packages("simsurv")}else{library(simsurv)}
  #library(simsurv)
  if (!require("MASS")){install.packages("MASS")}else{library(MASS)}
  #library(MASS)
  if (!require("deepTL")){devtools::install_github("SkadiEye/deepTL")}else{library(deepTL); print("SurvDeepFIT's dependencies successfully loaded. Congratulations!")}
  #library(deepTL)
  if (!require("survivalmodels")){install.packages("survivalmodels")}else{library(survivalmodels)}
  #library(survivalmodels)
  if (!require("randomForestSRC")){install.packages("randomForestSRC")}else{library(randomForestSRC)}
  #library(randomForestSRC)
  if (!require("survivalsvm")){install.packages("survivalsvm")}else{library(survivalsvm)}
  #library(survivalsvm)
  if (!require("reticulate")){install.packages("reticulate")}else{library(reticulate)}
  #library(reticulate)
  if (!require("Hmisc")){install.packages("Hmisc")}else{library(Hmisc)}
  #library(Hmisc)
  if (!require("xgboost")){install.packages("xgboost")}else{library(xgboost)}
  #library(xgboost)
  if (!require("pheatmap")){install.packages("pheatmap")}else{library(pheatmap)}
  if (!require("stringr")){install.packages("stringr")}else{library(stringr)}
  if (!require("timeROC")){install.packages("timeROC")}else{library(timeROC)}
  #library(pheatmap)
  conda_list = reticulate::conda_list()
  conda_path = conda_list[which(conda_list$name==use_condaenv),2]
  use_condaenv(conda_list[which(conda_list$name==use_condaenv),2])
  deep = deepsurv(formula = Surv(time,status)~sexF+age+trt,data = simsurvdata(50)) %>% try()
  if("try-error" %in% class(deep)) {
    if (str_detect(deep,"pycox")){
      print("Please install PyCox by survivalmodels::install_pycox() or reticulate::conda_install()")
      print("Suggested pycox version: 0.2.3")
    }else if (str_detect(deep,"pytorch")){
      print("Please install PyTorch by reticulate::py_install(envname  = use_condaenv,packages = \"pytorch==1.9.0\")")
      print("Suggested torch version: 1.9.0")
    }
  }else{
    print("DeepSurv and DeepHit successfully loaded. Congratulations!")
  }
}


#' Title Cox-PH log-likelihood calculation
#'
#' @param t_threshold The threshold time to calculate the risk set.
#' @param times Patients' event times.
#'
#' @return Returns the IDs of patients whose event times are longer than t_threshold.
#' @export
#'
risk.set <- function(t_threshold,times) {
  return(which(times >= t_threshold))
}



#' Title Calculate the Cox-PH Log-Likelihood Based on Event Status, Event Time, and Risk Prediction
#'
#' @param Status Patients' event statuses.
#' @param Times Patients' event times.
#' @param f_hat_y The product of X and Beta (X \%*\% Beta).
#'
#' @return Returns the log-likelihood of a Cox-PH model.
#' @export
#'
#'
loglik_coxph = function(Status,Times,f_hat_y){
  if (length(Status) != length(Times)){
    return("Check your input!")
  } else {
    loglik = 0
    Xbeta = f_hat_y %>% as.numeric() # 被替换成为 predict value的地方
    df= data.frame(status=Status,
                   times=Times,
                   yhat=f_hat_y)
    df$yhat[which(is.infinite(df$yhat))] = NA
    df$loglik = apply(df,1,function(x){
      status = x["status"]
      time = x["times"]
      riskset= risk.set(t_threshold = time,times = df$times)
      yhat_i = x["yhat"]
      yhat_riskset= df$yhat[riskset]
      return(as.numeric(status)*(yhat_i-log(sum(exp(yhat_riskset),na.rm = T))))
    })
    df$loglik[which(is.infinite(df$loglik))]=NA
    loglik = sum(df$loglik,na.rm = T)
  }
  if (is.infinite(loglik)){
    return(NA)
  }else{
    return(loglik)
  }
}



#' Title Internal Implementations of Machine Learning Models, Including RSF, Lasso-Cox, etc.
#'
#' @param method Name of the method:
#' - "random_forest" for Random Survival Forest;
#' - "survival_aft" for Accelerated Failure Time Model;
#' - "survival_cox" for Cox-PH Model;
#' - "DeepSurv" for DeepSurv;
#' - "DeepHit" for DeepHit;
#' - "Xgboost" for XGBoost;
#' - "Survival_SVM" for Survival Support Vector Machine;
#' - "lasso" for Lasso-Cox;
#' - "ensemble_dnnet" for SurvDeepFIT.
#'
#'
#' @param model.type Default is "survival" for survival data.
#' @param object A dnnetSurvInput object, created by deepTL::importDnnetSurv().
#' @param ... Other hyper-parameters for the machine learning models.
#' @author Shiyu Wan
#'
#' @return Trained machine learning models.
#' @export
#'
mod_permfit <- function(method, model.type, object, ...) {
  if (model.type == "survival"){
    if (method == "random_forest"){
      mod <- do.call(randomForestSRC::rfsrc,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "survival_aft"){
      df_train=cbind(object@y,object@e,object@x) %>% as.data.frame()
      colnames(df_train)=c("EventTime","EventStatus",
                           colnames(object@x))
      mod <- survival::survreg(Surv(EventTime,EventStatus)~.,data = df_train,dist = "weibull")
    } else if (method == "survival_cox") {
      df_train=cbind(object@y,object@e,object@x) %>% as.data.frame()
      colnames(df_train)=c("EventTime","EventStatus",
                           colnames(object@x))
      mod <- survival::coxph(Surv(EventTime,EventStatus)~.,data = df_train)
    } else if (method == "DeepSurv") {
      mod <- do.call(survivalmodels::deepsurv,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "DeepHit") {
      mod <- do.call(survivalmodels::deephit,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "Xgboost"){
      label_train = ifelse(object@e == 1,object@y,-object@y)
      df_train_xg = xgboost::xgb.DMatrix(data = as.matrix(object@x),label = label_train)
      xgb_args = list(...)
      xgb_params = xgb_args$params
      if (is.null(xgb_params)) {
        xgb_params = list()
      }
      if (!is.null(xgb_params$nrounds)) {
        if (is.null(xgb_args$nrounds)) {
          xgb_args$nrounds = xgb_params$nrounds
        }
        xgb_params$nrounds = NULL
      }
      if (is.null(xgb_args$nrounds)) {
        xgb_args$nrounds = 100L
      }
      if (!is.null(xgb_args$objective)) {
        xgb_args$objective = NULL
      }
      xgb_train_args = names(formals(xgboost::xgb.train))
      xgb_param_args = setdiff(names(xgb_args), xgb_train_args)
      for (xgb_arg in xgb_param_args) {
        xgb_params[[xgb_arg]] = xgb_args[[xgb_arg]]
        xgb_args[[xgb_arg]] = NULL
      }
      xgb_params$objective = "survival:cox"
      xgb_args$params = xgb_params
      mod =do.call(xgboost::xgb.train,
                   appendArg(xgb_args,"data",df_train_xg,TRUE))
    } else if (method == "Survival_SVM"){
      mod <- do.call(survivalsvm::survivalsvm,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "DnnSurv") {
      mod <- do.call(survivalmodels::dnnsurv,
                     appendArg(appendArg(list(...), "formula", Surv(y, e) ~ ., TRUE),
                               "data", data.frame(x = object@x, y = object@y, e = object@e), TRUE))
    } else if (method == "lasso") {
      lasso_family <- "cox"
      cv_lasso_mod <- glmnet::cv.glmnet(object@x, cbind(time = object@y,status = object@e), family = lasso_family)
      mod <- glmnet::glmnet(object@x, cbind(time = object@y,status = object@e), family = lasso_family,
                            lambda = cv_lasso_mod$lambda[which.min(cv_lasso_mod$cvm)])
    } else if (method == "ensemble_dnnet") {
      mod <- do.call(deepTL::ensemble_dnnet, appendArg(list(...), "object", object, TRUE))
    } else if (method == "dnnet") {
      spli_obj = splitDnnet(object, 0.8)
      mod <- do.call(dnnet, appendArg(appendArg(list(...), "train", spli_obj$train, TRUE),
                                      "validate", spli_obj$valid, TRUE))
    } else {
      return("Survival Analysis Not Applicable")
    }
  } else {
    if(method == "ensemble_dnnet") {
      mod <- do.call(ensemble_dnnet, appendArg(list(...), "object", object, TRUE))
    } else if (method == "random_forest") {
      mod <- do.call(randomForest::randomForest,
                     appendArg(appendArg(list(...), "x", object@x, TRUE), "y", object@y, TRUE))
    } else if (method == "lasso") {
      lasso_family <- ifelse(model.type == "regression", "gaussian",
                             ifelse(model.type == "binary-classification", "binomial", "cox"))
      cv_lasso_mod <- glmnet::cv.glmnet(object@x, object@y, family = lasso_family)
      mod <- glmnet::glmnet(object@x, object@y, family = lasso_family,
                            lambda = cv_lasso_mod$lambda[which.min(cv_lasso_mod$cvm)])
    } else if (method == "linear") {
      if(model.type == "regression") {
        mod <- stats::lm(y ~ ., data.frame(x = object@x, y = object@y))
      } else if(model.type == "binary-classification") {
        mod <- stats::glm(y ~ ., family = "binomial", data = data.frame(x = object@x, y = object@y))
      } else {
        mod <- survival::coxph(survival::Surv(y, e) ~ ., data = data.frame(x = object@x, y = object@y, e = object@e))
      }
    } else if (method == "svm") {
      if(model.type == "regression") {
        mod <- e1071::tune.svm(object@x, object@y, gamma = 10**(-(0:4)), cost = 10**(0:4/2),
                               tunecontrol = e1071::tune.control(cross = 5))
        mod <- mod$best.model
      } else if(model.type == "binary-classification") {
        mod <- e1071::tune.svm(object@x, object@y, gamma = 10**(-(0:4)), cost = 10**(0:4/2),
                               tunecontrol = e1071::tune.control(cross = 5))
        mod <- svm(object@x, object@y, gamma = mod$best.parameters$gamma, cost = mod$best.parameters$cost, probability = TRUE)
      } else {
        return("Not Applicable")
      }
    } else if (method == "dnnet") {
      spli_obj <- splitDnnet(object, 0.8)
      mod <- do.call(dnnet, appendArg(appendArg(list(...), "train", spli_obj$train, TRUE),
                                      "validate", spli_obj$valid, TRUE))
    } else {
      return("Not Applicable")
    }
  }
  return(mod)
}

#' Title Internal Implementations of Predicting Algorithms for Different Machine Learning Models
#'
#' @param mod Trained machine learning models, returned by mod_permfit().
#' @param object Test dataset, a dnnetSurvInput object created by deepTL::importDnnetSurv().
#' @param method Type of the machine learning models, e.g., "survival_cox" for Cox-PH model.
#' @param model.type Default is "survival" for survival data.
#' @author Shiyu Wan
#'
#' @return Risk predictions for patients in the test dataset.
#' @export
#'
predict_mod_permfit <- function(mod, object, method, model.type) {
  if(model.type == "regression") {
    if(!method %in% c("linear", "lasso")) {
      return(predict(mod, object@x))
    } else if(method == "linear") {
      return(predict(mod, object@x)[, "s0"])
    } else {
      return(predict(mod, object@x)[, "s0"])
    }
  } else if(model.type == "binary-classification") {
    if(method %in% c("dnnet")) {
      return(predict(mod, object@x)[, mod@label[1]])
    } else if(method == "ensemble_dnnet") {
      return(predict(mod, object@x)[, mod@model.list[[1]]@label[1]])
    } else if(method == "random_forest") {
      return(predict(mod, object@x, type = "prob")[, 1])
    } else if (method == "lasso") {
      return(1 - predict(mod, object@x, type = "response")[, "s0"])
    } else if (method == "linear") {
      return(1 - predict(mod, data.frame(x = object@x, y = object@y), type = "response"))
    } else if (method == "svm") {
      return(attr(predict(mod, object@x, decision.values = TRUE, probability = TRUE),
                  "probabilities")[, levels(object@y)[1]])
    }
  } else if(model.type == "survival") {
    if(method %in% c("ensemble_dnnet", "dnnet")) {
      pred = predict(mod, object@x)
      return(list(pred,pred))
    } else if(method == "random_forest") {
      pred = as.numeric(predict(mod, data.frame(x = object@x))$predicted)
      return(list(pred,pred))
    } else if (method == "lasso") {
      pred = predict(mod, object@x)[, "s0"] %>% as.numeric() %>% exp()
      return(list(pred,log(pred)))
    } else if (method == "survival_cox") {
      pred = predict(mod, data.frame(object@x), type = "risk")
      return(list(pred,log(pred)))
    } else if (method == "survival_aft") {
      pred = -as.numeric(predict(mod, data.frame(object@x), type = "lp"))
      return(list(pred,pred))
    } else if (method == "Survival_SVM") {
      pred = -as.numeric(predict(mod, data.frame(x = object@x))$predicted)
      return(list(pred,pred))
    } else if (method == "DeepSurv" | method == "DeepHit" | method == "DnnSurv") {
      pred = predict(mod, data.frame(x = object@x),type = "risk")
      return(list(pred,log(pred)))
    } else if (method == "Xgboost"){
      pred = predict(mod,newdata = object@x)
      return(list(pred,log(pred)))
    }
  } else {
    return("Not Applicable")
  }
}


#' Title Internal Method to Calculate Harrell's C-Index
#'
#' @param Status 	A numeric vector of patients' survival status: 1 = event, 0 = censored.
#' @param Times  A numeric vector of patients' survival times.
#' @param f_hat_y  A numeric vector of patients' survival risk predictions.
#'
#' @return A numeric value representing Harrell's C-Index.
#' @export
#'
Cindex = function(Status,Times,f_hat_y){
  Inf_fhaty = which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)
  if (sum(which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)) == 0){
    rcorr_cens = rcorr.cens(f_hat_y,Surv(Times,Status))[1]
    Cindex = 1 - as.numeric(rcorr_cens)
    return(Cindex)
  }else{
    rcorr_cens = rcorr.cens(f_hat_y[-Inf_fhaty],Surv(Times[-Inf_fhaty],Status[-Inf_fhaty]))[1]
    Cindex = 1 - as.numeric(rcorr_cens)
    return(Cindex)
  }
}



#' Title Internal Method to Calculate Time-Dependent AUC
#'
#' @param Status A numeric vector of patients' survival status: 1 = event, 0 = censored.
#' @param Times A numeric vector of patients' survival times.
#' @param f_hat_y A numeric vector of patients' survival risk predictions.
#'
#' @return The average of time-dependent AUCs at the 25\%, 50\%, and 75\% quantiles of the event time.
#' @export
#'
timedAUC = function(Status,Times,f_hat_y){
  Inf_fhaty = which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)
  if (sum(which(is.infinite(f_hat_y) == T | is.na(f_hat_y) == T)) == 0){
    death_case = which(Status == 1)
    choose_time = Times[death_case]
    time_range = quantile(choose_time,c(0.25,0.5,0.75))
    tdAUC = timeROC(T = Times, delta = Status,
                    marker = f_hat_y,weighting = "marginal",cause = 1,
                    times = time_range)
    return(mean(tdAUC$AUC))
  }else{
    new_time = Times[-Inf_fhaty]
    new_status = Status[-Inf_fhaty]
    new_f_hat_y = f_hat_y[-Inf_fhaty]
    death_case = which(new_status == 1)
    choose_time = new_time[death_case]
    time_range = quantile(choose_time,c(0.25,0.5,0.75))
    tdAUC = timeROC(T = new_time, delta = new_status,
                    marker = new_f_hat_y,weighting = "marginal",cause = 1,
                    times = time_range)
    return(mean(tdAUC$AUC))
  }
}

#' Title  Internal Method for Measuring Prediction Performance Differences in Survival Models
#'
#' @param model.type Default is "survival" for survival data.
#' @param y_hat A numeric vector of patients' survival risk predictions.
#' @param y_hat0 A numeric vector of patients' survival risk predictions based on permuted data (details explained in the manuscript).
#' @param object Validation dataset, a dnnetSurvInput object created by deepTL::importDnnetSurv().
#' @param y_max Inner upper boundary for y in binary classification; can be ignored when model.type is "survival".
#' @param y_min Inner lower boundary for y in binary classification; can be ignored when model.type is "survival".
#' @param y_hatcoxl A numeric vector of patients' survival risk predictions returned by Cox-related models, e.g., Cox, SurvDeepFIT, XGBoost.
#' @param y_hat0coxl A numeric vector of patients' survival risk predictions returned by Cox-related models based on permuted data.
#'
#' @return A numeric vector consisting of differences in C-index, Cox's partial log-likelihood, and average time-dependent AUC, based on observed and permuted data.
#' @export
#'
log_lik_diff <- function(model.type, y_hat, y_hat0, object,
                         y_max = 1-10**-10, y_min = 10**-10,
                         y_hatcoxl,y_hat0coxl) {
  if(model.type == "regression") {
    return((object@y - y_hat)**2 - (object@y - y_hat0)**2)
  } else if(model.type == "binary-classification") {
    y_hat <- ifelse(y_hat < y_min, y_min, ifelse(y_hat > y_max, y_max, y_hat))
    y_hat0 <- ifelse(y_hat0 < y_min, y_min, ifelse(y_hat0 > y_max, y_max, y_hat0))
    return(-(object@y == levels(object@y)[1])*log(y_hat) - (object@y != levels(object@y)[1])*log(1-y_hat) +
             (object@y == levels(object@y)[1])*log(y_hat0) + (object@y != levels(object@y)[1])*log(1-y_hat0))
  } else if(model.type == "survival") {
    return(c(Cindex(Status = object@e,Times = object@y,f_hat_y = y_hat)-
               Cindex(Status = object@e,Times = object@y,f_hat_y = y_hat0),
             loglik_coxph(Status = object@e,Times = object@y,f_hat_y = y_hatcoxl)-
               loglik_coxph(Status = object@e,Times = object@y,f_hat_y = y_hat0coxl),
             timedAUC(Status = object@e,Times = object@y,f_hat_y = y_hat)-
               timedAUC(Status = object@e,Times = object@y,f_hat_y = y_hat0)))
  } else {
    return("Not Applicable")
  }
}

##### Update 07/03/2024: To solve the identifiability issue, relative risk will be centered around zero.

#' Title Centralizing the Relative Risk Prediction of SurvDeepFIT Around 0
#'
#' @param mod A fitted dnnet or dnnetEnsemble object.
#' @param object The training set corresponding to the above model.
#'
#' @return A centralized model with the last bias updated, and the event times of the training set.
#' @export
#'
model.centralize = function(mod,object){
  if (class(mod)[1] == "dnnet"){
    Uncentralized_Risk = predict(mod,object@x)
    mod@bias[[length(mod@bias)]] = mod@bias[[length(mod@bias)]] - mean(Uncentralized_Risk,na.rm = T)
  }else if (class(mod)[1] == "dnnetEnsemble"){
    for (i in 1:length(mod@model.list)) {
      Uncentralized_Risk_i = predict(mod@model.list[[i]],object@x)
      mod@model.list[[i]]@bias[[length(mod@model.list[[i]]@bias)]] = mod@model.list[[i]]@bias[[length(mod@model.list[[i]]@bias)]] -mean(Uncentralized_Risk_i)
    }
  }
  return(list(mod,object))
}

#' Title Internal Method to Predict Patient's Survival Based on Relative Risk Prediction and Breslow's Estimator of Baseline Cumulative Hazard
#'
#' @param centralized_model The output of the function model.centralize().
#' @param test_object Test dataset, a dnnetSurvInput object.
#' @param time_point If specified, it corresponds to the specific time points where you want to predict the patient's survival; otherwise, the survival probability at each unique failure time in the training dataset will be provided.
#'
#' @return A data.frame with each patient's survival probability prediction.
#' @export
#'
predict_surv_df = function(centralized_model,test_object,time_point = NULL){
  mod = centralized_model[[1]]
  train_object = centralized_model[[2]]
  fail_time = train_object@y[which(train_object@e==1)]
  unique_time_in_train = unique(train_object@y)
  sort_time_in_train = sort(unique_time_in_train)
  unique_fail_time = unique(fail_time)
  sort_fail_time = sort(unique_fail_time)
  train_risk_score = predict(mod,train_object@x)
  test_risk_score = predict(mod,test_object@x)
  # Estimate deltaH0(T_i) based on ## Time-to-Event Prediction with Neural Networks and Cox Regression
  dH0 = c()
  for (i in 1:length(sort_time_in_train)) {
    T_i = sort_time_in_train[i]
    Death_at_T_i = which(train_object@y == T_i & train_object@e == 1)
    D_i = length(Death_at_T_i)
    Riskset_i = risk.set(t_threshold = T_i,times = train_object@y)
    #dH_0_Ti = D_i/sum(exp(train_risk_score)[Riskset_i])
    dH_0_Ti = D_i/sum(exp(train_risk_score)[which(train_object@y > T_i)])
    dH0 = append(dH0,dH_0_Ti)
  }
  if (is.null(time_point)){
    # Only Calculate H0(t) for all the failure time in the training set
    # If specific time is not provided, then all the failure time points in the training set will be considered
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train<=t))])
      H0 = append(H0,H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0*exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time,3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test,surv_df)
    }
    return(Surv_df_test)
  }else{
    sort_fail_time = sort(time_point)
    H0 = c()
    for (i in 1:length(sort_fail_time)) {
      t = sort_fail_time[i]
      H0_t = sum(dH0[1:last(which(sort_time_in_train<=t))])
      H0 = append(H0,H0_t)
    }
    Surv_df_test = data.frame()
    for (i in 1:length(test_risk_score)) {
      surv_df = t(as.data.frame(exp(-H0*exp(test_risk_score[i]))))
      colnames(surv_df) = round(sort_fail_time,3)
      rownames(surv_df) = i
      Surv_df_test = rbind(Surv_df_test,surv_df)
    }
    return(Surv_df_test)
  }
}

##### 4. Permfit Survival #####

#' Title PermFIT: A permutation-based feature importance test extended to survival analysis
#'
#' @param train A training dataset, which is a dnnetSurvInput object created by deepTL::importDnnetSurv().
#' @param validate A validation dataset is required when k_fold = 0.
#' @param k_fold K-fold cross-fitting. Default is k_fold = 5. If preferring not to use the cross-fitting strategy, set k_fold to zero.
#' @param n_perm Number of permutations repeated for each feature, which is R in the manuscript's algorithm 1. Default is n_perm = 100.
#' @param pathway_list For categorical variables, dummy variables should be created in advance. For variables with more than 2 categories, all the corresponding dummy variables should be included in the pathway_list and permuted simultaneously to calculate the permutation feature importance.
#' @param method Name of the method:
#' - "random_forest" for Random Survival Forest;
#' - "survival_aft" for Accelerated Failure Time Model;
#' - "survival_cox" for Cox-PH Model;
#' - "DeepSurv" for DeepSurv;
#' - "DeepHit" for DeepHit;
#' - "Xgboost" for XGBoost;
#' - "Survival_SVM" for Survival Support Vector Machine;
#' - "lasso" for Lasso-Cox;
#' - "ensemble_dnnet" for SurvDeepFIT.
#' @param shuffle If NULL, the data will be shuffled for cross-fitting; if random shuffle is not desired, provide a vector of numbers for cross-fitting indices.
#' @param ... Additional parameters passed to each method.
#'
#' @return A PermFIT object, with importance = permutation feature importance for all continuous and binary features, and imp_block = permutation feature importance for all categorical features with more than two groups.
#' @export
#'
permfit_survival <- function(train, validate = NULL, k_fold = 5,
                             n_perm = 100, pathway_list = list(),
                             method = c("ensemble_dnnet", "random_forest",
                                        "lasso", "linear", "svm", "dnnet",
                                        "survival_cox","survival_aft","DeepSurv",
                                        "DeepHit","Xgboost","Survival_SVM",
                                        "DnnSurv")[8],
                             shuffle = NULL,
                             ...){
  n_pathway <- length(pathway_list)
  # Number of observations
  n <- dim(train@x)[1]
  # Number of Covariates
  p <- dim(train@x)[2]

  # Decide model type based on input data
  if(class(train) == "dnnetSurvInput") {
    model.type <- "survival"
  } else if(class(train) == "dnnetInput") {
    if(is.factor(train@y)) {
      model.type <- "binary-classification"
    } else {
      model.type <- "regression"
    }
  } else {
    stop("'train' has to be either a dnnetInput or dnnetSurvInput object.")
  }

  # No Cross validation, need a validation set
  if(k_fold == 0) {
    if(is.null(validate)){
      stop("A validation set is required when k = 0. ")
    }
    # Sample Size of Validation set
    n_valid <- dim(validate@x)[1]
    # fit model
    mod = mod_permfit(method = method,model.type = model.type,object = train,...)
    # Risk Score
    f_hat_x <- predict_mod_permfit(mod = mod,object = validate,
                                   method = method,model.type = model.type)
    valid_ind <- list(1:length(validate@y))
    y_pred <- f_hat_x

    # n_pathway: for categorical variables，several dummy variable will be created
    # need to permutation all such dummies simultaneously


    if(n_pathway >= 1) {
      p_score <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorel <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorea <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      for(i in 1:n_pathway) {
        for(l in 1:n_perm){
          # Validate Dataset
          set.seed(i*1000+l)
          x_i <- validate@x
          # 只permutation一个Block
          x_i[, pathway_list[[i]]] <- x_i[, pathway_list[[i]]][sample(n_valid), ]
          Perm_Validate = importDnnetSurv(x = x_i,
                                          y = validate@y,
                                          e = validate@e)
          pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                       object = Perm_Validate,method = method)
          diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                            y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                            y_hat0coxl = pred_i[[2]],object = train_spl$train)
          p_score[l, , i] = diff[1]
          p_scorel[l, , i] = diff[2]
          p_scorea[l, , i] = diff[3]
        }
      }
    }

    # PERMFIT:Core
    p_score2 <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2l <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2a <- array(NA, dim = c(n_perm, n_valid, p))

    for (i in 1:p) {
      for (l in 1:n_perm) {
        set.seed(i*1000+l)
        x_i <- validate@x
        x_i[, i] <- x_i[, i][sample(dim(x_i)[1])]
        Perm_Validate = importDnnetSurv(x = x_i,
                                        y = validate@y,
                                        e = validate@e)
        pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                     object = Perm_Validate,method = method)
        diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                          y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                          y_hat0coxl = pred_i[[2]],object = train_spl$train)
        p_score2[l, , i] = diff[1]
        p_score2l[l, , i] = diff[2]
        p_score2a[l, , i] = diff[3]
      }
    }
  } else { # K-Fold
    valid_ind <- list()
    # Shuffle the sample
    if(is.null(shuffle)) {
      shuffle <- sample(n)
    }
    n_valid <- n
    y_pred <- numeric(length(train@y))

    ## n_pathway：for dummies
    if(n_pathway >= 1){
      p_score <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorel <- array(NA, dim = c(n_perm, n_valid, n_pathway))
      p_scorea <- array(NA, dim = c(n_perm, n_valid, n_pathway))
    }

    # initialize permutation score matrix
    p_score2 <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2l <- array(NA, dim = c(n_perm, n_valid, p))
    p_score2a <- array(NA, dim = c(n_perm, n_valid, p))

    valid_error <- numeric(k_fold)

    # K-Fold:: Start
    for(k in 1:k_fold){
      # Splid the data
      train_spl <- splitDnnet(train, shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)])
      # 记录所用的验证集的观测编号
      valid_ind[[k]] <- shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)]
      # Typo: Valid is the 4 of 5 folds
      mod = mod_permfit(method = method,model.type = model.type,object = train_spl$valid,...)
      #print(mod)
      # Risk Score
      f_hat_x <- predict_mod_permfit(mod = mod,object = train_spl$train,
                                     method = method,model.type = model.type)

      #print(f_hat_x)
      # valid_error[k] <- sum(log_lik_diff(model.type, f_hat_x, f_hat_x, train_spl$train))

      y_pred[valid_ind[[k]]] <- f_hat_x[1]

      if(k == 1) {
        final_model <- mod
      } else if(method == "ensemble_dnnet") {
        final_model@model.list <- c(final_model@model.list, mod@model.list)
        final_model@loss <- c(final_model@loss, mod@loss)
        final_model@keep <- c(final_model@keep, mod@keep)
      }

      # To Be Done
      if(n_pathway >= 1){
        for(i in 1:n_pathway) {
          for(l in 1:n_perm) {
            x_i <- train_spl$train@x
            x_i[, pathway_list[[i]]] <- x_i[, pathway_list[[i]]][sample(dim(x_i)[1]), ]
            Perm_Validate = importDnnetSurv(x = x_i,
                                            y = train_spl$train@y,
                                            e = train_spl$train@e)
            pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                         object = Perm_Validate,method = method)
            diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                              y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                              y_hat0coxl = pred_i[[2]],object = train_spl$train)
            p_score[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[1]
            p_scorel[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[2]
            p_scorea[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[3]
          }
        }
      }
      # Permfit::Core
      # p: p the covariate

      for (i in 1:p) {
        for (l in 1:n_perm) {
          set.seed(i*1000+l)
          x_i <- train_spl$train@x
          x_i[, i] <- x_i[, i][sample(dim(x_i)[1])]
          Perm_Validate = importDnnetSurv(x = x_i,
                                          y = train_spl$train@y,
                                          e = train_spl$train@e)
          pred_i = predict_mod_permfit(mod = mod,model.type = model.type,
                                       object = Perm_Validate,method = method)
          diff=log_lik_diff(model.type = model.type ,y_hat = f_hat_x[[1]],
                            y_hat0 = pred_i[[1]],y_hatcoxl = f_hat_x[[2]],
                            y_hat0coxl = pred_i[[2]],object = train_spl$train)
          p_score2[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[1]
          p_score2l[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[2]
          p_score2a[l,shuffle[floor((k-1)*n/k_fold+1):floor(k*n/k_fold)], i] = diff[3]
        }
      }
      print(k)
    }
    mod <- final_model
    # valid_error <- sum(valid_error)/n_valid
  }

  # Importance Score: C-Index
  if(is.null(colnames(train@x))) {
    imp <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    imp <- data.frame(var_name = colnames(train@x))
  }
  imp$importance <- apply(apply(p_score2, 2:3, mean), 2, mean, na.rm = TRUE)
  imp$importance_sd <- sqrt(apply(apply(p_score2, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
  imp$importance_pval <- 1 - stats::pnorm(imp$importance/imp$importance_sd)
  if(n_perm > 1) {
    imp$importance_sd_x <- apply(apply(p_score2, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    imp$importance_pval_x <- 1 - stats::pnorm(imp$importance/imp$importance_sd_x)
  }

  # Importance Score: likelihood
  if(is.null(colnames(train@x))) {
    impl <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    impl <- data.frame(var_name = colnames(train@x))
  }
  impl$importance <- apply(apply(p_score2l, 2:3, mean), 2, mean, na.rm = TRUE)
  impl$importance_sd <- sqrt(apply(apply(p_score2l, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
  impl$importance_pval <- 1 - stats::pnorm(impl$importance/impl$importance_sd)
  if(n_perm > 1) {
    impl$importance_sd_x <- apply(apply(p_score2l, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    impl$importance_pval_x <- 1 - stats::pnorm(impl$importance/impl$importance_sd_x)
  }

  # Importance Score: tdAUC
  if(is.null(colnames(train@x))) {
    impa <- data.frame(var_name = paste0("V", 1:p))
  } else  {
    impa <- data.frame(var_name = colnames(train@x))
  }
  impa$importance <- apply(apply(p_score2a, 2:3, mean), 2, mean, na.rm = TRUE)
  impa$importance_sd <- sqrt(apply(apply(p_score2a, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
  impa$importance_pval <- 1 - stats::pnorm(impa$importance/impa$importance_sd)
  if(n_perm > 1) {
    impa$importance_sd_x <- apply(apply(p_score2a, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
    impa$importance_pval_x <- 1 - stats::pnorm(impa$importance/impa$importance_sd_x)
  }

  imp_block <- data.frame()
  if(n_pathway >= 1) {
    if(is.null(names(pathway_list))) {
      imp_block <- data.frame(block = paste0("P", 1:n_pathway))
    } else {
      imp_block <- data.frame(block = names(pathway_list))
    }
    imp_block$importance <- apply(apply(p_score, 2:3, mean), 2, mean, na.rm = TRUE)
    imp_block$importance_sd <- sqrt(apply(apply(p_score, 2:3, mean), 2, stats::var, na.rm = TRUE)/n_valid)
    imp_block$importance_pval <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd)
    if(n_perm > 1) {
      imp_block$importance_sd_x <- apply(apply(p_score, c(1, 3), mean), 2, stats::sd, na.rm = TRUE)
      imp_block$importance_pval_x <- 1 - stats::pnorm(imp_block$importance/imp_block$importance_sd_x)
    }
  }

  return(new("PermFIT", model = list(impl,impa), importance = imp, block_importance = imp_block,
             validation_index = valid_ind, y_hat = y_pred))
}
