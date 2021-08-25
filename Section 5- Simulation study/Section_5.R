rm(list=ls())
library("survival")
source("checks.R")
source("functions.R") ## for sampling, etc.
require(plyr)
require(dplyr)
options(digits=10)
# set.seed(123)
# N = 50

#' Examples
#' rtest(get_par(N=100), vis=T)   -- standard test with correct specification
#' rtest(get_par(K1=0, K2=0, gamma=2), misspecify = T, vis = T) --- misspecified model

fit_cox <- function(adoption) {
  # Fit a Cox model.
  adoption_cox = adoption
  if("cox_data" %in% names(adoption)) {
    adoption_cox = adoption$cox_data
  }
  sur = with(adoption_cox, Surv(start, stop, event))
  X = adoption_cox$X
  fit = coxph(formula = sur ~ X)               
  beta_hat = coef(fit)          
  ps_hat =  exp(X * beta_hat) / sum (exp(X * beta_hat))
  return(list(beta_hat=beta_hat, ps_hat=ps_hat))
}


# Various Test statistics:
#  DID = difference-in-differences
#  DR = Doubly robust version
Tstat_DID_DR <- function(treated_i, adoption, Y) {
  par = adoption$par
  
  X = adoption$X
  n = length(X)
  t1 = adoption$t1
  pp = break_pre_post(t1, par)
  pre = c(pp$pre)
  post = c(pp$post)
  rplus = mean(par$rho^post/ (1-par$rho))
  rminus = mean(par$rho^pre / (1-par$rho))
  
  Sn = Tstat_DID(treated_i, adoption, Y)
  
  Sn - (n / (n-1)) * par$gamma * (rplus - rminus) * (mean(X) - X[treated_i])
}

# Fast implementation
Tstat_DID_DR_fast <- function(treated_i, adoption, Y) {
  #' 1. estimate gamma
  #'
  par_hat = fast_hats(adoption, Y)
  rho_hat = par_hat$rho
  
  gamma_hat = par_hat$gamma
  gamma_hat = adoption$par$gamma
  
  t1 = adoption$t1
  pp = break_pre_post(t1, adoption$par)
  pre = c(pp$pre)
  post = c(pp$post)
  
  rplus = mean(rho_hat^post/ (1-rho_hat))
  rminus = mean(rho_hat^pre / (1-rho_hat))
  
  Sn = Tstat_DID(treated_i, adoption, Y)
  n = nrow(Y)
  X = adoption$X
  stopifnot(n == length(X))
  Sn - (n / (n-1)) * gamma_hat * (rplus - rminus) * (mean(X) - X[treated_i])
}


Tstat_DID <- function(treated_i, adoption, Y) {
  t1 = adoption$t1
  pp = break_pre_post(t1, adoption$par)
  pre = c(pp$pre)
  post = c(pp$post)
  Y11 = Y[treated_i, post]; Y10 = Y[treated_i, pre]
  Y01 = colMeans(matrix(Y[-treated_i, post], ncol=length(post))); 
  Y00 = colMeans(matrix(Y[-treated_i, pre], ncol=length(pre)))
  (mean(Y11 - Y01) -  mean(Y10 - Y00))
}

# Matching test statistic (not used in main paper.)
Tstat_MATCH = function(treated_i, adoption, Y) {
  t1 = adoption$t1
  pp = break_pre_post(t1, adoption$par)
  post = c(pp$post)
  
  X = adoption$X
  J = order(abs(X - X[treated_i]))[2]
  mean(Y[treated_i, post] - Y[J, post])
}

#' SCM = Synthetic Control Test statistic (Eq. (17) in main paper)
Tstat_SCM <- function(t1, unit1, Y, X){
  # create matrices from panel data that provide inputs for synth(
  t1_index = to_index(t1)
  df = data.frame(y=c(), unit=c(), time=c())
  n = length(X)
  for(i in 1:n) {
    M = matrix(c(rep(i, nT), Y[i, ], 1:nT, rep(X[i], nT), c(0, head(Y[i, ], nT-1))), ncol=5, byrow=F)
    colnames(M) <- c("unit", "y", "time", "x", "yprev")
    df = rbind(df, M)  
  }
  
  post = seq(t1_index + 1, nT)
  pre = seq(1, t1_index)
  dataprep.out<-
    dataprep(
      foo = df,
      predictors = c("x", "yprev"),
      predictors.op = "mean",
      dependent = "y",
      unit.variable = "unit",
      time.variable = "time",
      treatment.identifier = unit1,
      controls.identifier = setdiff(c(1:n), unit1),
      time.predictors.prior = pre,
      time.optimize.ssr = pre,
      # unit.names.variable = c("state")
      time.plot = 1:nT
    )
  
  synth.out <- synth(dataprep.out, verbose = F)
  w = as.numeric(synth.out$solution.w)
  
  y1pre_synth = as.numeric(t(Y[-unit1, pre]) %*% w)
  y1post_synth = as.numeric(t(Y[-unit1, post]) %*% w)
  
  pre_gap <-  as.numeric(Y[unit1, pre]) - y1pre_synth
  post_gap = as.numeric(Y[unit1, post]) - y1post_synth
  
  mean(post_gap^2) / mean(pre_gap^2)
}

Tstat_mle <- function(adoption, Y) {
   out = CI_mle(t1, treated_i, Y, X)
   out$point
}

# ALL_Tstat = list("DID"=Tstat_DID, "DID_DR"=Tstat_DID_DR,  "DID_DR_fast"=Tstat_DID_DR_fast)
# Could also include matching statistic: "MATCH"=Tstat_MATCH)

#' Main function to implement a randomization test.
#' @param par Parameter values for the simulation. See functions.R and get_par().
#' @param misspecify Whether to use a misspecified model or not (see 5.2 Robustness section)
#' @param Tstat_list List of test statistic functions to use. Paper uses only DID.
#' @param nreps No. of resamples in the randomization test (more is better.)
#' 
rtest <- function(par, tau=0, nreps=1000, 
                  Tstat_list=ALL_Tstat,
                  misspecify=FALSE, vis=FALSE) {
  #' 
  t0 = proc.time()[3]
  Results = matrix(0, nrow=0, ncol=4)
  colnames(Results)  = c("uniform", "phat", "oracle", "method")
  summarize_results = function(RR) {
    df = as.data.frame(RR)
    S = ddply(df, .(method), summarize, 
              unif=round(100*mean(uniform), 3), 
              phat=round(100*mean(phat), 3), 
              oracle=round(100*mean(oracle), 3),
              num=length(uniform))
    S$method = names(Tstat_list)[S$method]
    return(S)
  }
  
  n = par$N # number of units
  nT = par$nT # number of time points
  if(misspecify) {
    print(">>> Using misspecified adoption model to sample treatment times.")
  }

  for(j in 1:length(Tstat_list)) {
    Tstat = Tstat_list[[j]]
    
    for(irep in 1:nreps) {
      # 1. Sample data
      X = sample_X(par)
      adopt = sample_Adoption(X, par)
      if(misspecify) {
        adopt = sample_misspecified_Adoption(X, par)
        par$Tmax = adopt$Tmax  # need to adjust because Tmax is calculated with baseline model.
        par$TIME = seq(0, par$Tmax, length.out=par$nT)
        adopt$par = par
      }
      
      I = adopt$I # first adopter unit
      t1 = adopt$t1 # T_(1), time of first adoption
      Y = sample_Y(adopt, par, tau=tau)
      
      # 2. Compute all test stat values (for speed)
      tvals = sapply(1:n, function(i) Tstat(i, adopt, Y))
      tobs = tvals[I]
      
      # 3. Propensity scores.
      true_ps =  adopt$ps
      unif_ps = rep(1/n, n)
      
      fit = tryCatch({ 
        fit_cox(adopt$cox_data)
      }, 
      error=function(cond) {
        return(list(ps_hat=unif_ps))
      }, 
      warning=function(cond) {
        return(list(ps_hat=unif_ps))
      })
      
      # 2. Perform test.
      ps_hat = fit$ps_hat / sum(fit$ps_hat)
      
      
      if(any(is.na(ps_hat))) {
        ps_hat = unif_ps # if problem with estimation, fall back to permutation test.
      }
      
      t1 = proc.time()[3]
      Results = rbind(Results, 
                      c(two_sided_test(tobs, tvals, unif_ps, alpha=0.05),
                        two_sided_test(tobs, tvals, ps_hat, alpha=0.05),
                        two_sided_test(tobs, tvals, true_ps, alpha=0.05),
                        j))
      
      if(t1 - t0 > 15) {
        par(mfrow=c(1, 1))
        print(sprintf("rtest(): Iter = %d/%d -- Misspecified Model? %s", irep, nreps, misspecify))
        S = summarize_results(Results)
        print(S)
        if(vis) {
          par(mfrow=c(2, 1))
          plot(X, ps_hat, xlab="X", 
               ylab="estimated ps", pch=20, col="blue")
          abline(0, 1, col="red")
          plot(X, true_ps, pch=20)
          abline(0, 1, col="red")
        }
        t0 = t1
      } 
    }
  }
  
  S = summarize_results(Results)
  print(S)
  return(S)
}

get_sim_params = function(sim_id, make_it_quick) {
  if(sim_id==1) {
    A = expand.grid(gamma=c(0, 0.5, 1, 2, 5), 
                    N = c(25, 50, 100), 
                    tau = c(0))
    # For quick replication.
    if(make_it_quick) {
      A = expand.grid(gamma=c(0, 0.5, 2), 
                      N = c(25), 
                      tau = c(0))
    }
    return(A)  
  } else if(sim_id==2) {
    B = expand.grid(gamma=c(0, 0.5, 1, 2, 5), 
                    N = c(25, 50, 100), 
                    tau = c(0.25, 0.5))
    # For quick replication.
    if(make_it_quick) {
      B = expand.grid(gamma=c(0, 0.5, 2), 
                      N = c(25), 
                      tau = c(0.5))
    }
    
    return(B)  
  } else if(sim_id==3 || sim_id==4) {
    
    C = expand.grid(gamma=c(0, 2, 5),
                    K1 = c(0, 1, 2),
                    K2 = c(0, 1, 2),
                    N = c(50), 
                    tau=c(0))
    
    # For quick replication.
    if(make_it_quick) {
      C = expand.grid(gamma=c(0, 5),
                      K1 = c(0, 2),
                      K2 = c(0, 2),
                      N = c(25),
                      tau=c(0))
    }
  } else {
    stop("Not correct Simulation ID.")
  }
}

#' Function to reproduce Tables from the paper.
#' 
Reproduce_Table = function(table_id, nreps=1000, make_it_quick=TRUE) {
  
  A = NULL # simulation settings
  cols = c("unif", "phat", "oracle")
  
  ALL_Tstat = list("DID"=Tstat_DID)
  if(table_id==4) {
    ALL_Tstat = list("DID_DR"=Tstat_DID_DR_fast) # doubly robust 
  }
  misspecify = FALSE
  if(table_id >=3 ) {
    misspecify=TRUE
  }
  if(make_it_quick) {
    nreps = ifelse(nreps > 500, 500, nreps)
  }
  
  Results = matrix(0, nrow=0, ncol=length(cols) + 3)
  colnames(Results) = c("γ", "n","τ", cols)
  if(table_id >= 3) {
    colnames(Results) = c("γ", "k1","k2", cols)
  }
  
  A = get_sim_params(table_id, make_it_quick)
  print(paste(">> Generating table ", table_id, " with ", nreps, " randomization reps. Quick Mode? ", make_it_quick))
  
  for(i in 1:nrow(A)) {
    print(paste(">.. Running simulation setting ", i, " / ", nrow(A), "...<"))
    tau = A$tau[i];
    N = A$N[i]
    g = A$gamma[i]
   

    par=NULL
    k1=NULL; k2=NULL
    if(table_id < 3) {
      par = get_par(N=N, gamma=g)
      print(paste(">>> γ = ", g, " (model param.), n = ", N, "(#units)",", τ = ", tau, "(treatment effect)"))
      
    } else {
      par = get_par(N=N, gamma=g, K1 = A$K1[i], K2 = A$K2[i])
      k1 = par$K1
      k2 = par$K2
      print(paste(">>> γ = ", g, " (model param.), n = ", N, "(#units)",",k1, k2 = ", k1,k2))
    }
    
    out = rtest(par=par, nreps = nreps, Tstat_list = ALL_Tstat, misspecify = misspecify, tau = tau, vis=F)
    M = as.numeric(as.vector((as.matrix(out)[, -c(1, 5)])))
    if(table_id < 3) {
      Results = rbind(Results, c(g, N, tau, M))
    } else {
      Results = rbind(Results, c(g, k1, k2, M))
    }
  }
  
  print(paste("% rejection rates over ", nreps, "repetitions."))
  
  return(Results)
  
}
