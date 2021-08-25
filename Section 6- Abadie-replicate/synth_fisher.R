rm(list=ls())
library(plyr)
library(survival)
require(xtable)

## Loads adoption_DataList
## Contains several specifications of the DATA.
source("synth_fisher_DATA.R")


load("smoking_with_dems.rdata")
smoking$unit.num = NULL
cols = c(3:8)  # interesting covariates.

#' Used to de-trend X_it
time_trends = function() {
  ## Covariates
  X = matrix(0, nrow=0, ncol=35)
  before1980 = seq(1, 10)
  X = rbind(X, ddply(smoking, .(year), summarize, mean(cigsale))[-before1980, 2]) # baseline is 1980
  X = rbind(X, ddply(smoking, .(year), summarize, mean(lnincome))[-before1980, 2])
  X = rbind(X, ddply(smoking, .(year), summarize, mean(age15to24))[-before1980, 2])
  X = rbind(X, ddply(smoking, .(year), summarize, mean(retprice))[-before1980, 2])
  X = rbind(X, ddply(smoking, .(year), summarize, mean(unemploy))[-before1980, 2])
  X = rbind(X, ddply(smoking, .(year), summarize, mean(dems))[-before1980, 2])
  rownames(X) = colnames(smoking)[cols]
  B_trend = matrix(0, nrow=0, ncol=1)
  for(j in 1:nrow(X)) {
    Xj = X[j, ]
    z = 1:ncol(X)
    B_trend = rbind(B_trend, coef(lm(Xj ~ z))[2])
  }
  rownames(B_trend) = rownames(X)
  
  #' B_trend = 
  #' cigsale   -2.3869101487
  #' lnincome   0.0309519193
  #' age15to24 -0.0008508221
  #' retprice  15.3226861961
  #' unemploy  -0.0337556561
  #' dems      -0.5573361715
  return(B_trend)
}

#' Main function.
#' adoption = data frame with treatment adoption data (state, when)
#' @return Matrix with propensity scores (P(I_1|..)) for every state.
#' 
synth_fisher = function(adoption, vars=c("lnincome", "retprice"), only.CA=FALSE, verbose=FALSE) {
  
  # fits Xt  ~ t to take out time-effect.
  X_trend = time_trends()
  
  # Week format
  week = sapply(adoption$when, function(s) {
    yr = as.numeric(strsplit(as.character(s), "/")[[1]][2])
    mo = as.numeric(strsplit(as.character(s), "/")[[1]][1])
    12*(yr-1980) + mo   # baseline is 1980.
  })
  
  #' Define covariates.
  #' 
  #' X_ij = covariate j of unit i at time of treatment.
  #' X_ij is defined in terms of "1980 values" where we difference out common time trends.
  #' 
  X = matrix(0, nrow=0, ncol=6)
  colnames(X) = colnames(smoking)[cols]
  stopifnot(all(vars %in% colnames(X)))
  
  # Adjust Xit for time trends.
  for(i in 1:length(AllStates)) {
    yr  = as.numeric(strsplit(as.character(adoption[i,]$when), "/")[[1]][2])
    st = as.character(adoption[i, ]$state)
    # risk_set = AllStates[which(week >= week[i])]
    
    # st_data = subset(smoking, state==st & year==yr)[, c(3:7, 9)]
    x_it = as.numeric(subset(smoking, state==st & year==yr)[, cols])
    x_it = x_it - as.numeric(X_trend) * (yr - 1980) 
    
    X = rbind(X, x_it)
  }
  rownames(X) = NULL
  # Update adoption data.
  adoption = cbind(adoption, X)
  head(adoption)
  # state    when  cigsale  lnincome age15to24    retprice unemploy     dems
  # Alabama 05/2004 145.1858  9.765749 0.1755933 -44.3444687 6.510136 84.80464
  
  # Change from (01/1990) -> week format
  adoption$when = week
  head(adoption)
  #' state    when  cigsale  lnincome  age15to24     retprice  unemploy     dems
  #' Alabama  293  145.1858  9.765749  0.1755933  -44.3444687  6.510136   84.8046
  #' ...
  
  status = rep(1, length(AllStates))
  status[which(AllStates=="Missouri")] = 0
  adoption$event = status
  
  surv = with(adoption, Surv(when, event))
  f = as.formula(paste("surv ~ ", paste(vars, collapse="+")))
  
  out = coxph(f, data=adoption)
  if(verbose) {
    print(sprintf("## Model ##"))
    print(sprintf("AIC = %.2f", AIC(out)))
    print(summary(out))
    print("## ##")
  }
  #
  var_ids = as.numeric(sapply(vars, function(v) which(colnames(adoption)==v)))
  
  X = as.matrix(adoption[, var_ids])
  stopifnot(all(names(coef(out)) == colnames(X)))
  
  # hats.
  yhat = exp(X %*% as.numeric(coef(out)))
  ps_hat = yhat / sum(yhat)
  rownames(ps_hat) = AllStates
  if(only.CA) {
    i = which(rownames(ps_hat)=="California")
    return(c(ps_hat[i, ], AIC(out)))
  }
  
  ord = rev(order(ps_hat))
  M = data.frame(a = rep(0, 13))
  # matrix(0, nrow=13, ncol=6)
  for(j in 1:3) {
    j1 = 13 * (j-1) + 1
    j_index = ord[seq(j1, j1 + 12)]
    ps_j = round(as.numeric(ps_hat[j_index]), 4)
    names_j = rownames(ps_hat)[j_index]
    
    M = cbind(M, data.frame(State=names_j, PS=ps_j))
  }
  
  M$a = NULL
  rownames(M) = NULL
  
  as.matrix(M)
}

#' Try all possible combinations of models with X1, ... X6
#' 
#' @return Kx3 matrix that contains (pvalue, AIC, #vars) at each row
#' 
single_DATA_analysis = function(adoption, verbose=FALSE) {
  # out = synth_fisher(adoption, vars = c("lnincome", "retprice"), FALSE)

  ## All models
  print("- Checking all models with 1-6 variables...")
  pvalues = matrix(0, nrow=0, ncol=3)
  colnames(pvalues) = c("pvalue", "AIC", "vars")
  
  for(num_var in 1:6) {
    # All num_var models
    models = t(combn(colnames(smoking)[cols], num_var))
    for(i in 1:nrow(models)) {
      m = models[i, ]
      # print(sprintf("Checking model"))
      # print(m)
      out = synth_fisher(adoption, vars = m, only.CA = TRUE)
      pvalues = rbind(pvalues, c(out, num_var))
    }
  }
  
  rownames(pvalues) = NULL
  #' pvalues = MATRIX (K x 3)
  #'    pvalue  AIC  #vars.
  #'         ....
  return(as.data.frame(pvalues))
}

# single_DATA_analysis(adoption_Data)

#' Analyzes all data adoption specifications in "adoption_DataList
#' @return Results matrix
#'       (dataset_id, %non_reject, mean_AIC, bestYES)
#' "bestYes focuses on Q1 of models in terms of AIC (i.e., best 25% of models)
#' 
full_DATA_analysis = function() {
  Results = matrix(0, nrow=0, ncol=4)
  colnames(Results) = c("dataSpec", "nonReject", "meanAIC", "isBest")
  K = length(adoption_DataList)
  
  for(i in 1:K) {
    print(sprintf("Testing data specification %d / %d ", i, K))
    dat = adoption_DataList[[i]]
    pvals = single_DATA_analysis(dat)
    nonRej = mean(pvals[, 1] >= 0.05)
    avgAIC = mean(pvals[, 2])
    M = pvals
    M = cbind(rep(i, nrow(pvals)), M)
    Results = rbind(Results, M)
    
    # best = which(pvals[, 2] <= quantile(pvals[, 2], .25))
    # nonRej_best = mean(pvals[best, 1] >= 0.05)
    # avgAIC_best = mean(pvals[best, 2])
    # Results = rbind(Results, c(i, nonRej_best, avgAIC_best, 1))
    
    print(Results)
    save(Results, file="SynthFisher_Results.rda")
  }
  return(Results)
}


paper_analysis = function() {
  # single_DATA_analysis(adoption_Data_2, colnames(smoking)[cols])
  all_xnames = colnames(smoking)[cols]
  out = synth_fisher(adoption_Data, vars=c("lnincome", "retprice"),verbose=T)
  out
  xtable(out, include.rownames=FALSE)
  
  out = synth_fisher(adoption_Data_2, vars=c("lnincome", "retprice"), verbose=T)
  out
  xtable(out, include.rownames=FALSE)
  
  ## Single data analysis
  pvals = single_DATA_analysis(adoption_Data)
  pvals[which.min(pvals$AIC), ]
}

paper_appendix = function() {
  #
  load("SynthFisher_Results.rda")  
  Results = as.data.frame(Results)
  all = Results
  head(all)
  
  require(ggplot2)
  g = ggplot(data=all, aes(x=pvalue, y=AIC))
  g = g + geom_point() + xlab("pvalue") + ylab("average AIC")
  plot(g)
}
