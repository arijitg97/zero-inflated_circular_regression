library(circular)
library(CircStats)
library(splines2)
library(numDeriv)
library(parallel)


# l2 norm
l2 = function(v){
  sqrt(sum(v^2))
}

#projected data
x_pr = function(beta, x){
  as.vector(x %*% (beta/l2(beta)))
}

beta = function(phi) c(1/sqrt(1 + t(phi) %*% phi), as.numeric(1/sqrt(1 + t(phi) %*% phi))*phi)

jac.beta = function(phi){
  a1 = -as.numeric(1/(1 + t(phi) %*% phi)^(3/2))*phi
  a2 = as.numeric(1/sqrt(1 + t(phi)%*%phi))*diag(1, 2) - as.numeric(1/(1 + t(phi) %*% phi)^(3/2))*(phi %*% t(phi))
  return(matrix(rbind(a1, a2), 3, 2, byrow = T))
}


cdist = function(a, b){
  delta = abs(a - b)     
  pmin(delta, 2*pi - delta)
}


# Example data (angles in radians)
d = read.csv("dat.csv")  # Reading "dat.csv" after saving it in the current working directory

x0 = cbind(d[, 2], sin(d[, 3]), cos(d[, 3])) # Covariates
x = apply(x0, 2, function(x) x - mean(x))    # Centred covariates
y = d[, 1]                                   # Response


################################################



int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}


E.vmsp = function(x, y, phi, gamma, mu, df = d.f){
  b = int_b(x_pr(beta(phi), x), df = df)
  delta = y - mu - 2*atan(as.vector(b %*% gamma))
  return(-sum(cos(delta)))
}


lc.vmsp = function(x, y, init, df = d.f, tol = 1e-3, max.iter = 200){
  theta.0 = init
  diff = 1
  iter = 0
  q = ncol(x) - 1
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[3:(2+q)]), x), df = df)
    
    gamma_new = optim(theta.0[-c(1:(2+q))], function(gamma) E.vmsp(x, y, theta.0[3:(2+q)], gamma, theta.0[1], df = df), method = "BFGS", control = list(maxit = 2000))$par
    phi_new = optim(theta.0[3:(2+q)], function(phi) E.vmsp(x, y, phi, gamma_new, theta.0[1], df = df), method = "BFGS", control = list(maxit = 2000))$par
    
    b = int_b(x_pr(beta(phi_new), x), df = df)
    del = y - 2*atan(as.vector(b %*% gamma_new))
    C = mean(cos(del))
    D = mean(sin(del))
    
    mu_new = atan2(D, C) 
    kappa_new = A1inv(sqrt(C^2 + D^2))
    theta.0 = c(mu_new, kappa_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)     
  }
}

loglik.vmsp = function(theta, x, y, df = d.f){
  q = ncol(x) - 1
  b = int_b(x_pr(beta(theta[3:(2+q)]), x), df = df)
  mu = theta[1] + 2*atan(as.vector(b %*% theta[-c(1:(2+q))]))
  sum(log(dvm(y, mu, theta[2])))
}



vmsp = function(x, y, n.knots = 3, m = 200){
  #Fit Fisher-Lee model
  fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type = 'c-l'))
  
  d.f = n.knots + 2
  # Initial values
  beta0 = fl$coefficients
  phi0 = beta0[-1]/beta0[1]
  phi.init = matrix(phi0, nrow = m, ncol = length(phi0), byrow = TRUE)
  gamma.init = matrix(rnorm(m*d.f, 0, 0.3), nrow = m, ncol = d.f, byrow = TRUE)
  init = cbind(as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)
  
  # Parallel computation 
  ncores = detectCores()
  cl = makeCluster(ncores - 1)
  clusterEvalQ(cl, {library(splines2)
    library(CircStats)
    library(circular)
  })
  clusterExport(cl, c("l2", "beta", "x_pr", "int_b", "E.vmsp", "lc.vmsp", "loglik.vmsp", "x", "y", "init", "d.f"),
                envir = environment())
  res = parLapply(cl, 1:m, function(i) {
    tryCatch({
      theta.m = lc.vmsp(x, y, init[i,], df = d.f)
      if (is.null(theta.m)) stop("did not converge")
      ll = loglik.vmsp(theta.m, x, y, df = d.f)
      c(theta.m, ll)
    }, error = function(e) NULL)
  })
  res = Filter(Negate(is.null), res)
  res = if (length(res) > 0) do.call(rbind, res) else NULL
  stopCluster(cl)
  
  ll = res[, ncol(res)]
  theta.m = res[, -ncol(res)]
  theta.hat = theta.m[which.max(ll),]
  aic = 2*length(theta.hat) - 2*loglik.vmsp(theta.hat, x, y, df = d.f)
  n = length(y)
  bic = length(theta.hat)*log(n) - 2*loglik.vmsp(theta.hat, x, y, df = d.f)
  q = ncol(x) - 1
  est = c(theta.hat[1:2], beta(theta.hat[3:(2+q)]), theta.hat[-c(1:(2+q))])
  cov = solve(hessian(function(theta) -loglik.vmsp(theta, x, y, df = d.f), theta.hat))
  se = sqrt(diag(cov))
  cov.beta = jac.beta(theta.hat[3:(2+q)]) %*% cov[3:(2+q), 3:(2+q)] %*% t(jac.beta(theta.hat[3:(2+q)]))
  se.beta = sqrt(diag(cov.beta))
  se_est = c(se[1:2], se.beta, se[-c(1:(2+q))])
  b = est[3:(3+q)]
  g = est[-c(1:(3+q))]
  
  # Estimated mean direction
  B_0 = int_b(x_pr(b, x), df = d.f)
  f_0 = as.vector(B_0 %*% g)
  mu.vmsp = (est[1] + 2*atan(f_0)) %% (2*pi)
  
  # Sum of residuals (to be computed for non-zero responses and corresponding estimated mean directions)
  id = (y != 0)
  e1 = sum(cdist(y[id], mu.vmsp[id]))
  
  #CRPS
  n.sim = 1e6
  set.seed(1)
  sim = runif(n.sim, 0, 2*pi)
  sim1 = runif(n.sim, 0, 2*pi)
  sim2 = runif(n.sim, 0, 2*pi)
  crp11 = sapply(1:n, function(i){
    2*pi * mean(cdist(sim, y[i]) * (CircStats::dvm(sim, mu.vmsp[i], est[2])))
  })
  crp12 = sapply(1:n, function(i){
    4*pi^2 * mean(cdist(sim1, sim2) *
                    (CircStats::dvm(sim1, mu.vmsp[i], est[2])) *
                    (CircStats::dvm(sim2, mu.vmsp[i], est[2])))
  })
  crps1 = crp11 - 0.5*crp12
  
  result = list(est[1:(3+q)], se_est[1:(3+q)], aic, bic, e1, mean(crps1))
  names(result)=c("Estimate", "Std.Error", "AIC", "BIC", "Sum_Residual", "Avg.CRPS")
  names(result[[1]]) = names(result[[2]]) = c("mu", "kappa", paste0("beta_", 1:(q+1)))
  result = lapply(result, round, 4)
  result$theta = theta.hat
  result$cov = cov
  result$d.f = d.f
  return(result)
}

t1 = Sys.time()
fit = vmsp(x,y)
t2 = Sys.time()
t2 - t1
fit


####### LOOCV #########


int1_b = function(x.new, x.ref, df = d.f, degree = 2){
  ibs_ref = ibs(x.ref, df = df, degree = degree)
  ibs_new = predict(ibs_ref, x.new)
  ibs_0   = predict(ibs_ref, 0)
  sweep(ibs_new, 2, ibs_0, "-")
}



n.sim = 1e5
set.seed(1)
sim = runif(n.sim, 0, 2*pi)
sim1 = runif(n.sim, 0, 2*pi)
sim2 = runif(n.sim, 0, 2*pi)

n = length(y)
cv_vmsp = numeric(n)
q = ncol(x) - 1

t1 = Sys.time()
suppressWarnings({
  for(i in 1:n){
  
  test.idx = i
  train.idx = setdiff(1:n, i)
  
  # Split predictors and response
  x.train = x[train.idx, , drop = FALSE]
  x.test  = x[test.idx, , drop = FALSE]
  
  y.train = y[train.idx]
  y.test  = y[test.idx]
 
  fit = vmsp(x.train, y.train)
  theta = fit$theta
  b = beta(theta[3:(2+q)])  
  g = theta[-c(1:(2+q))]
  
  xp.train = x_pr(b, x.train)
  xp.test  = x_pr(b, x.test)
  
  xp.test = pmin(pmax(xp.test, min(xp.train)),max(xp.train))
  
  if(length(unique(round(xp.train, 8))) <= 3){
    pred = NA
  } else {
    B_0 = try(int1_b(xp.test, xp.train, df = fit$d.f))
    if(inherits(B_0, "try-error")){
      pred = NA
    } else {
      f_0 = sum(g*B_0[1, ])
      pred = (theta[1] + 2*atan(f_0)) %% (2*pi)
    }
}
  if(all(is.na(pred))){
    cv_vmsp[i] = NA
  } else {
    crp1 = 2*pi * mean(cdist(sim, y.test)*(CircStats::dvm(sim, pred, theta[2])))
    crp2 = 4*pi^2 * mean(cdist(sim1, sim2)*(CircStats::dvm(sim1, pred, theta[2]))*(CircStats::dvm(sim2, pred, theta[2])))
    crp = crp1 - 0.5*crp2
    cv_vmsp[i] = mean(crp)
   }
  }
})
t2 = Sys.time()
t2 - t1

# Average CV error
round(mean(cv_vmsp), 4)






















