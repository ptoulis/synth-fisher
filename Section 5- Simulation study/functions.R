expit = function(u) {
  z = (u < -50) * (-50) + (u > 50) * 50 + (u > -50 & u < 50) * u
  exp(z) / (1+exp(z))
}

shape_misspecification = function() {
  A = expand.grid(k1=c(0, 2), k2=c(0, 2))
  par(mfrow=c(2, 2))
  par(mar=rep(4, 4))
  N = 250
  X = sample_X(list(N=N))
  
  for(i in 1:nrow(A)) {
    a = as.numeric(A[i,])
    par = get_par(N=N, K1=a[1], K2=a[2])
    n = length(X)
    m =  (1 - expit(2*a[1]*(X +8))) +  expit(2*a[2]*(X - 8))
    
    mu_X  = par$intercept  - 0.25*m
    sample_T = function() exp(mu_X + 0.4 * rnorm(n))
    
    Ti = sample_T()

    cnt = rep(0, n)
    for(k in 1:50000) { 
      # j = which.min(mu_X - rexp(n)) 
      j = which.min(sample_T())
      cnt[j] = cnt[j] + 1
    }
    true_ps = cnt / sum(cnt)
    J = order(X)
    plot(X[J], true_ps[J], pch=20, cex=.5,
         main=sprintf("k1= %.1f, k2=%.1f", a[1], a[2]),
         ylab="early adopter probability", xlab="X")
  }
}

get_par = function(N=50, nT=1000, intercept=1, beta=0.5, beta2=1, K1=2,K2=2,
                   rho=0.2, delta=0.5, gamma=0.5, sigma=0.15) {

  # For X-sampling.
  #' Figure out what Tmax should be and how to discretize given 
  #' the specified parameter values (rate_i = exp(icpt + beta *))
  #' 
  adopt_times = replicate(10000, { 
    Ti = sample_Adoption(sample_X(list(N=N)), 
                         list(beta=beta, intercept=intercept), 
                         retTi = TRUE) 
    as.numeric(quantile(Ti, probs = c(0.15)))
  })
  Tmax = 1*max(adopt_times)
  TIME = seq(0, Tmax, length.out=nT)
  
  # For Y sampling.
  A  = diag(nT-1)
  A = cbind(0, A)
  A = rbind(A, 0)
  C = solve(diag(nT) - rho * A)
  W_t =  matrix(sqrt(seq(1, nT)), nrow=N, ncol=nT, byrow=T)  # 1, ...sqrt(t)
  
  W = delta * W_t # design matrix with sqrt(t_i) values. To speed up sample_Y
  
  # matrices
  # print(paste(">>> N=", N, "Tmax=",round(Tmax, 2), 
            #  "beta=",beta, "beta2=", beta2, "K1=", K1, " K2=", K2, " gamma=",gamma, "delta=",delta, "sigma=",sigma))
  # print(paste(">> Done.."))
  
  list(N=N, nT=nT, TIME=TIME, Tmax=Tmax, 
       beta=beta, intercept=intercept,beta2=beta2,K1=K1,K2=K2,
       gamma=gamma, delta=delta, sigma=sigma, rho=rho, 
       W=W, W_t=W_t, C=C)  
}

#' Break timepoint into "pre" and "post" period.
break_pre_post = function(t1, par) {
  t1_index = max(which(par$TIME <= t1))
  stopifnot(t1_index < par$nT)
  pre = seq(1,  t1_index)
  post = seq(t1_index + 1, par$nT)
  return(list(pre=pre, post=post))
}

#' Sample covariates. Nx1 vector.
sample_X = function(par) {
  10*runif(par$N, min=-1, max=1)
  # sort(sample(c(-1, 0, 1), size=N, replace=T, prob=c(0.7, 0.2, 0.1)))
  # sample(c(-1, 0, 1), size=N, replace=T, prob=c(0.7, 0.2, 0.1))
}

#' Sample adoption based on Cox model.
#' 
sample_Adoption <- function(X, par, retTi = FALSE) {
  # Sample event times from Cox model.
  beta = par$beta # true parameter for Cox model
  kIntercept = par$intercept  # increase this to decrease Ti

  eta = exp(kIntercept + X * beta)
  true_ps = eta / sum(eta)
  n = length(X); stopifnot(n==par$N)
  Ti = rexp(n, rate=eta) # all event times
  
  if(retTi) return(Ti)
  
  Tmax = par$Tmax
  events = (Ti <= Tmax)  # events (TRUE=event happened)
  time0 = rep(0, n)
  time1 = events * Ti  + (1-events) * Tmax # treatment time
  
  nT = par$nT
  W_x = matrix(X, nrow=n, ncol=nT, byrow=F) # X matrix over time. 
  
  #y = sample_Y(min(Ti), which(Ti == min(Ti)), tau, X, delta=0.5, gamma=0.5)
  df = list(cox_data=data.frame(X=X, start=time0, 
                                stop=time1, event=events, 
                                ps=true_ps), 
            t1=min(time1), 
            I=which.min(time1), 
            X=X, ps=true_ps,
            W_x=W_x, 
            par=par)
  return(df)
}

old__calculate_misspecified_PS = function() {
  par = get_par(N=1000)
  n = length(X)
  mu_X = par$intercept + (X+10)^2 * par$beta2
  sample_T = function() exp(mu_X + 0.5 * rnorm(n))
  
  cnt = rep(0, length(X))
  for(k in 1:100000) { 
    # j = which.min(mu_X - rexp(n)) 
    j = which.min(sample_T())
    cnt[j] = cnt[j] + 1
  }
  ps = cnt/sum(cnt)
  par(mfrow=c(1, 1))
  plot(X, ps, pch=20)
  
  
  obj = function(theta, retvals=F) {
    A = theta[1] 
    B = theta[2]
    C = theta[3]
    D = theta[4]
    pshat = (X < C) * (A + B * (X + 10)) + (X>C) * D
    if(retvals) return(pshat)
    log(sum(ps-pshat)^2)
  }
  
  out = optim(par=rep(0, 4), obj, method="BFGS")
  out
  pshat = obj(out$par, retvals=T)
  pshat = obj(c(0.015, -0.013/2, -7, 0), T)
  plot(X, ps, col="red", pch=20)
  points(X, pshat, col="blue", pch=20, cex=0.2)
}

sample_misspecified_Adoption <- function(X, par, retTi = FALSE) {
  #' Uses accelerated failure model (AFT) instead of Cox.
  #' 
  stopifnot(c("K1", "K2") %in% names(par))
  beta = par$beta # true parameter for Cox model
  kIntercept = par$intercept  # increase this to decrease Ti
  n = length(X)
  
  # mu_X = kIntercept + (X>0) * (X-1)^2 * par$beta2 + (X<0) * (X+6)^2 * par$beta2
  # mu_X = kIntercept + (X-5)^2 * par$beta2 - (abs(X+5) < 1) * 8
  # mu_X = kIntercept + (X > 0) * (X-5)^2 * par$beta2  + (X < 0) * (X+10)^2 * par$beta2  
  # mu_X = kIntercept + abs(cos(2 * pi * (X - mean(X))/sd(X))) * par$beta2
  # Ti = exp(mu_X - rexp(n))  # Accelerated failure time model.
  # mu_X  = kIntercept  + (X+10)^2 * par$beta2 # - X * 0.001 * expit(5*X) * par$beta2
  quants = c(-8, 8)
  # as.numeric(quantile(X, probs=c(0.1, 0.9)))
  

  a_K =  (1 - expit(2*par$K1*(X - quants[1]))) +  expit(2*par$K2*(X - quants[2]))
  mu_X  = kIntercept  - a_K * par$beta2# - X * 0.001 * expit(5*X) * par$beta2
  sample_T = function() exp(mu_X + 0.4 * rnorm(n))
  
  Ti = sample_T()
  if(retTi) return(Ti)
  # Calculate true PS through simulation.
  cnt = rep(0, n)
  for(k in 1:10000) { 
    # j = which.min(mu_X - rexp(n)) 
    j = which.min(sample_T())
    cnt[j] = cnt[j] + 1
  }
  true_ps = cnt / sum(cnt)
  #true_ps = rep(1/n, n)
  # true_ps = sapply(1:n, function(i) {
  #   d = (mu_X[i] - mu_X[-i])
  #   mean(sapply(U_Exp, function(j) prod(pexp(j-d))))
  # })
  # true_ps =  true_ps/sum(true_ps)
  
  # Tmax = par$Tmax
  Tmax = quantile(Ti, 0.2)
  
  nT = par$nT
  events = (Ti <= Tmax)  # events (TRUE=event happened)
  time0 = rep(0, n)
  time1 = events * Ti  + (1-events) * Tmax # treatment time
  
  W_x = matrix(X, nrow=n, ncol=nT, byrow=F) # X matrix.
  
  df = list(cox_data=data.frame(X=X, start=time0, 
                                stop=time1, event=events, 
                                ps=true_ps), 
            t1=min(time1), 
            I=which.min(time1), 
            X=X, ps=true_ps,
            beta=beta, 
            Tmax=Tmax,
            W_x=W_x, 
            par=par)
  return(df)
}


sample_Y = function(adoption, par, tau) {
  t1 = adoption$t1
  treated_i = adoption$I
  # unload params.
  nT = par$nT; rho = par$rho; delta=par$delta; gamma=par$gamma; sigma=par$sigma
  
  pp = break_pre_post(t1, par)
  post = pp$post

  n = par$N
  # nx T matrices of covariates
  W_tau = matrix(0, nrow=n, ncol=nT)
  W_tau[treated_i, post] = 1
  W_x = adoption$W_x
  
  # Errors.
  Eps = matrix(rnorm(n * nT, sd=sigma), nrow=n, ncol=nT)
  
  # Control potential outcomes. AR(1) model
  # Yit = rho * Yit_1 + gamma X + delta sqrt(t) + e_it
  Y0 = (par$W  + gamma * W_x + Eps) %*% par$C
  # Treated Outcomes.
  Y = Y0 + tau * W_tau
  
  return(Y)
}


one_sided_test = function(tobs, tvals, ps, alpha=0.05, tol=1e-15, ret_pval=F) {
  n = length(ps)
  stopifnot(length(tvals) == length(ps))
  stopifnot(length(tobs)==1)
  if(abs(sum(ps) - 1) > tol) {
    print(ps)
    print(paste("PS does not sum to one...", sum(ps)))
    stop("")
  }
  
  if(ret_pval) {
    pval = sum(ps * (tvals >= tobs))
    hist(tvals, breaks=20)
    abline(v=tobs, col="red", lwd=2)
    return(pval)
  }
  
  ord = order(tvals)
  A = cumsum(ps[ord])
  k = min(which(A >= (1-alpha)))  # t_{(k)} = c_n^ in Eq. (9)
  Tq = tvals[ord[k]]  ## t_{(k)}
  
  p1 = sum(ps[tvals > Tq])
  p0 = sum(ps[abs(tvals - Tq) < tol])
  p3 = ((alpha-p1)/p0)
  
  if(tobs > Tq) { 
    return(1) 
  } else if( tobs == Tq) { 
    # return((alpha-p1))/p0
    # return(p3)
    return(as.numeric(runif(1) <= p3))
  } else {
    return(0)    
  }
  
}

two_sided_test = function(tobs, tvals, ps, alpha=0.05, tol=1e-15) {
  n = length(ps)
  stopifnot(length(tvals) == length(ps))
  stopifnot(length(tobs)==1)
  if(abs(sum(ps) - 1) > tol) {
    print(ps)
    print(sum(ps))
    stop("")
  }
  
  ord = order(tvals)
  A = cumsum(ps[ord])
  k1 = min(which(A >= (1-alpha/2))) 
  k0 = min(which(A >= (alpha/2)))
  
  Tq1 = tvals[ord[k1]]  ## t_{(k)}
  Tq0 = tvals[ord[k0]]  ## 
  
  p1 = sum(ps[tvals > Tq1])
  p0 = sum(ps[tvals < Tq0])
  
  p1c = sum(ps[abs(tvals - Tq1) < tol])
  p0c = sum(ps[abs(tvals - Tq0) < tol])
  
  p31 = ((alpha/2 - p1)/p1c)
  p30 = ((alpha/2 - p0)/p0c)
  
  
  if(tobs > Tq1 | tobs < Tq0) { 
    return(1) 
  } else if( abs(tobs - Tq1) < tol) { 
    return(as.numeric(runif(1) <= p31))
    #
  } else if(abs(tobs - Tq0) < tol) { 
    return(as.numeric(runif(1) <= p30))
  } else {
    return(0)    
  }
  
}


mle = function(adoption, Y) {
  
  par = adoption$par
  t1 = adoption$t1
  X = adoption$X
  I = adoption$I
  N = par$N
  nT = par$nT
  
  pp = break_pre_post(t1, par)
  post = pp$post
  
  delta_coeff =  matrix(sqrt(seq(1, nT)), nrow=N, ncol=nT, byrow=T)
  #
  D1 = delta_coeff[, -1]
  D1 = as.vector(D1)
  gamma_coeff = matrix(X, nrow=N, ncol=nT, byrow=F)
  G1 = gamma_coeff[, -1]
  G1 = as.vector(G1)
  #
  #
  
  tau_coeff = matrix(0, nrow=N, ncol=nT)
  tau_coeff[I, post] = 1
  T1 = as.vector(tau_coeff[, -1])
  fit_0 = lm(as.vector(Y[,-1]) ~ 0 + as.vector(Y[,seq(1, nT-1)]) + 
               D1 + G1 + T1)
  
  se = summary(fit_0)$coefficients[4, 2]
  
  loglik = function(tau_par, ret_all=F) {
    T2 = matrix(0, nrow=N, ncol=nT)
    T2[I, post] = tau_par
    ynew = Y - T2
    
    # ynew = r * ynew_ + d * Delta g * Gamma + e
    ymin = as.vector(ynew[, seq(1, nT-1)])
    y = as.vector(ynew[, -1])
    
    fit = lm(y ~ 0+ ymin + D1 + G1)
    s = sd(fit$residuals)
    
    if(ret_all) {return(fit)}
    sum(dnorm(fit$residuals, sd=s, log=T))
  }
  
  out = optim(c(coef(fit_0)[4]), loglik, method="BFGS", control=list(fnscale=-1, maxit=1000))
  
  out2 = as.numeric(coef(loglik(out$par, ret_all=T)))
  
  return(list(rho=out2[1], delta=out2[2], gamma=out2[3], tau=as.numeric(out$par)))
}

fast_hats = function(adoption, Y) {
  # estimation of gamma
  par = adoption$par
  X  = adoption$X
  n = nrow(Y)
  stopifnot(n == length(X))
  # A0 = solve(diag(n) + matrix(1, nrow=n, ncol=n))
  A0 = diag(n) - (1 / (n+1)) * matrix(1, nrow=n, ncol=n)
  y = A0 %*% (Y[,1] - Y[1,1])
  x = A0 %*% (X - X[1])
  gamma_hat = coef(lm(y ~ x + 0))[1]
  
  #' #' 2. estimate tau
  #' t1 = adoption$t1
  I = adoption$I
  #' pp = break_pre_post(t1, par)
  #' post = c(pp$post)
  #' X = adoption$X
  #' J = order(abs(X - X[I]))[2]
  #' tau_hat = mean(Y[I, post] - Y[J, post])
  #' 
  #' 3. Estimate rho
  n1 = n - 1
  # A = solve(diag(n1) + matrix(1, nrow=n1, ncol=n1))
  A = diag(n1) - (1/(n1 + 1)) * matrix(1, nrow=n1, ncol=n1)
  
  nT = par$nT
  Ynew = A %*% Y[-I,2:nT]
  Y_1 = A %*%  Y[-I, seq(1, nT-1)]
  
  X1 = A %*% matrix(X[-I]-X[1], nrow=n1, ncol=nT-1)
  y = as.vector(Ynew-gamma_hat * X1)
   x=   as.vector(Y_1)
  rho_hat = coef(lm(y ~ x + 0))[1]
  
  return(list( rho=rho_hat, gamma=gamma_hat))
}

