
#' Tests whether sample_Adoption() has the desired freq properties.
#' 
CHECK_sample_Adoption = function(par=get_par(), nreps=1e4) {
  X = sample_X(par)
  d0 = sample_Adoption(X, par)
  ps = d0$ps  # theoretical propensity to be first adopter.
  
  t0 = proc.time()[3]
  cnt = rep(0, par$N)  ## frequency of who is the first adopter
  print("Comparing target propensity score (Ps) and empirical propensities..")
  for(j in 1:nreps) {
    d = sample_Adoption(X, par)  # sample adoption times conditional on X
    stopifnot(all(d$ps == ps)) # PS scores need to be fixed.
    I = d$I
    t1 = proc.time()[3]
    cnt[I] = cnt[I] + 1
    
    if(t1 - t0 > 5) {
      ps_sim = as.numeric(cnt / j)  # empirical propensities.
      mse = sqrt(mean((ps_sim - ps)^2))
      print(sprintf("i=%d / %d mse=%.5f", j, nreps, mse))
      par(mfrow=c(1, 1))
      plot(ps, ps_sim, pch=20)
      abline(0, 1, col="red")
      t0 = t1
    }
  }
  
  par(mfrow=c(1, 1))
  ps_sim = as.numeric(cnt/nreps)
  qqplot(ps_sim, ps, pch=20, ylab="target propensity score", xlab="empirical propensity score")
  abline(0, 1, col="red", lwd=2)
}

CHECK_sample_misspecified_Adoption = function(par=get_par(), nreps=1e4) {
  X = sample_X(par)
  d0 = sample_misspecified_Adoption(X, par)
  ps = d0$ps  # theoretical propensity to be first adopter.
  
  t0 = proc.time()[3]
  cnt = rep(0, par$N)  ## frequency of who is the first adopter
  print("Comparing target propensity score (Ps) and empirical propensities..")
  for(j in 1:nreps) {
    d = sample_misspecified_Adoption(X, par)  # sample adoption times conditional on X
    # stopifnot(all(d$ps == ps)) # PS scores need to be fixed.
    I = d$I
    t1 = proc.time()[3]
    cnt[I] = cnt[I] + 1
    
    if(t1 - t0 > 5) {
      ps_sim = as.numeric(cnt/j)  # empirical propensities.
      mse = sqrt(mean((ps_sim - ps)^2))
      print(sprintf("i=%d / %d mse=%.5f", j, nreps, mse))
      par(mfrow=c(1, 1))
      plot(ps, ps_sim, pch=20)
      abline(0, 1, col="red")
      t0 = t1
    }
  }
  
  par(mfrow=c(1, 1))
  ps_sim = as.numeric(table(out) / nreps)
  qqplot(ps_sim, ps, pch=20, ylab="target propensity score", xlab="empirical propensity score")
  abline(0, 1, col="red", lwd=2)
}


CHECK_adoption_times = function(par = get_par()) {
  ts = replicate(1000, {
    ad = sample_Adoption(sample_X(par), par)
    ad$t1
  })
  
  par(mfrow=c(1, 1))
  hist(ts, breaks=40)
  print(paste("Tmax = ", par$Tmax))
  y = sapply(ts, function(t1) max(break_pre_post(t1, par)$pre))
  print("Adoption time index. Should not be =1 frequently.")
  print(table(y))
}

#' Tests whether fit_cox() is estimating survival model well.
#' 
CHECK_fit_cox <- function(par=get_par(), verbose=T) {
  X = sample_X(par)
  d = sample_Adoption(X, par)
  colMeans(d$cox_data)
  fit = fit_cox(d$cox_data)
  mse =  sqrt(mean((fit$ps_hat - d$ps)^2))
  if(!verbose) return(mse)
  
  print(paste("MSE(ps - ps^) = ", mse))
  par(mfrow=c(1, 2))
  plot(d$ps, fit$ps_hat, pch=20)
  abline(0, 1, col="red")
  print(sprintf("beta=%.2f beta_hat=%.2f", par$beta, fit$beta_hat))
  plot(X, d$ps)
}

CHECK_misspecified_fit_cox <- function(par=get_par(), verbose=T) {
  X = sample_X(par)
  d = sample_misspecified_Adoption(X, par)
  
  fit = fit_cox(d$cox_data)
  mse =  sqrt(mean((fit$ps_hat - d$ps)^2))
  # if(!verbose) return(mse)
  
  print(paste("MSE(ps - ps^) = ", mse))
  par(mfrow=c(1, 2))
  plot(d$ps, fit$ps_hat, pch=20)
  abline(0, 1, col="red")
  print(sprintf("beta=%.2f beta_hat=%.2f", d$beta, fit$beta_hat))
  plot(X, d$ps)
}

CHECK_sample_Y = function(par=get_par()) {
  t0 = proc.time()[3]
  X = sample_X(par)
  ad = sample_Adoption(X, par)
  
  tau_val = -1
  Y = sample_Y(ad, par, tau=tau_val)
  J = order(abs(X - X[ad$I]))[2]
  
  par(mfrow=c(1, 1))
  plot(Y[ad$I, ], type="l", col="blue")
  lines(Y[J, ], col="red")
  print(paste("tau = ",   tau_val, " naive matching estimator = ", Tstat_MATCH(ad$I, ad, Y)))
  
  nreps=100
  y = replicate(nreps, { 
    #  set.seed(123)
    y1 = sample_Y(ad, par, 0)
    #   set.seed(123)
    # y2 = sample_Y_v2(ad, par0, 0)
    1
  })
  print(paste("sample_Y takes ", round((proc.time()[3] - t0)/ nreps, 3)," secs / call."))
}

#' Check whether the one-sided test has correct coverage.
#' 
CHECK_one_sided_test = function() {
  
  n = 50
  rej = replicate(200000, {
    u = runif(n/5)
    tvals = sample(u, size=n, replace=T)
    ps = abs(rt(n, df=2)); ps = ps/sum(ps)
    tobs = sample(tvals, size=1, prob = ps)
    
   #  ps_hat = 1/ps
    # ps_hat = ps_hat / sum(ps_hat)
    ps_hat = ps
    one_sided_test(tobs, tvals, ps_hat, alpha=0.05)
  })
  print(round(100*mean(rej), 2))
}


CHECK_two_sided_test = function() {
  n = 50
  
  rej = replicate(1e5, {
    
    u = runif(n/5)
    tvals = sample(u, size=n, replace=T)
    ps = abs(rt(n, df=2)); ps = ps/sum(ps)
    
    tobs = sample(tvals, size=1, prob = ps)
    two_sided_test(tobs, tvals, ps, alpha=0.05)
  })
  print(mean(rej))
}

CHECK_mle = function(tau, verb=F) {
  X = sample_X()
  dat = sample_Adoption(X)
  I = which.min(dat$stop)
  t1 = min(dat$stop); 
  Y = sample_Y(t1, I,  tau = tau, X = X, 
               rho=kRho, delta=kDelta, gamma=kGamma, sigma=1)

  time_0 = proc.time()[3]
  out = CI_mle(t1, I, Y, X)
  time_1 = proc.time()[3]
  if(verb) { print(paste("T1 - T0  = ", round(time_1 - time_0, 4), " secs"))
    print(out)
  }
  return(as.numeric(time_1 - time_0))
}


CHECK_DID = function(par=get_par()) {
  #' get_par(sigma=0, gamma=0) could be used to check true tau.
  #' 
  ok <<- c()
  
  e = replicate(1000, { 
    X = sample_X(par)
    dat = sample_Adoption(X, par)
    I = dat$I
    t1 = dat$t1
    n = length(X)
    
    tau = 1
    # t1  = first treatment time, I = first treated id.
    # Y = sample_Y(dat, par0, tau=tau)
    Y = sample_Y(dat, par, tau)
    
    pp = break_pre_post(t1, par)
    pre = pp$pre; post = pp$post
    stopifnot(length(post) > 0)
    
    ok <<- c(ok, length(pre) > 3)
    rplus = mean(par$rho^post/ (1-par$rho))
    rminus = mean(par$rho^pre / (1-par$rho))
    
    Sn = Tstat_DID(I, dat, Y)
    
    A = tau +  (n / (n-1)) * par$gamma * (rplus - rminus) * (mean(X) - X[I])
    Sn - A
  })

  #' e should be normal(0, sigma^2).
  #' 
  possible_s = exp(seq(log(sd(e)) -10, log(sd(e)) + 10, length.out = 1000))
  d0 = rnorm(5000)
  pvals = sapply(possible_s, function(s) {
    ks.test(d0, e/s)$p.value
  })
  print(paste("p-value of normality KS test: ", 
              round(max(pvals), 3), " for sigma = ", round(possible_s[which.max(pvals)], 3)))
  
  print(paste("mean = ", mean(e), " sd(e) = ", sd(e)))
  print(summary(e))
  par(mfrow=c(1, 1))
  hist(e, breaks=50)
  
  
  print(paste("% ok ", round(100*mean(ok), 2)))
}


CHECK_mle = function(nreps=10) {
  ##
  theta = runif(4, min=-1, max=1)
  par = get_par(rho=theta[1], delta=theta[2], gamma=theta[3])
  hats = matrix(0, nrow=0, ncol=4)
  colnames(hats) = c("rho", "delta", "gamma", "tau")
  times = c()
  for(i in 1:nreps) {
    X = sample_X(par)
    ad =sample_Adoption(X, par)
    Y = sample_Y(ad, par, tau=theta[4])  
    t0 = proc.time()[3]
    out = mle(ad, Y)
    t1 = proc.time()[3]
    times =c(times, t1-t0)
    hats = rbind(hats, unlist(out))
  }
  
  print(paste("average time for MLE computation ", mean(times), "sec."))
  par(mfrow=c(2, 2))
  for(j in 1:4) {
    hist(hats[,j])
    abline(v=theta[j], col="red", lwd=3)
  }
}


CHECK_fast_hats = function(nreps=10) {
  ##
  true_vals= matrix(0, nrow=0, ncol=2)
  
  times = c()
  hats = matrix(0, nrow=0, ncol=2)
  colnames(hats) = c("rho", "gamma")
  
  for(i in 1:nreps) {
    
    print(paste("Simulation - ", i, " / ", nreps, " .. "))
    theta = runif(3, min=-1, max=1)
    par = get_par(rho=theta[1], gamma=theta[2])
    
    X = sample_X(par)
    ad =sample_Adoption(X, par)
    Y = sample_Y(ad, par, tau=theta[3])  
    t0 = proc.time()[3]
    out = fast_hats(ad, Y)
    t1 = proc.time()[3]
    times =c(times, t1-t0)
    hats = rbind(hats, unlist(out))
    true_vals = rbind(true_vals, theta[1:2])
  }
  
  print(paste("average time for MLE computation ", mean(times), "sec."))
  par(mfrow=c(1, 2))
  for(j in 1:2) {
    plot(hats[,j], true_vals[,j], main=colnames(hats)[j])
    abline(0, 1, col="red")
  }
}