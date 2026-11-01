library(circular)
library(CircStats)
library(splines2)
library(numDeriv)
library(parallel)


# l2 norm
l2 = function(v){
  sqrt(sum(v^2))
}

# Projected data
x_pr = function(beta, x){
  as.vector(x %*% (beta/l2(beta)))
}

# Reparameterisation of beta
beta = function(phi) c(1/sqrt(1+t(phi)%*%phi), as.numeric(1/sqrt(1+t(phi) %*% phi))*phi)

# Jacobian
jac.beta = function(phi){
  a1 = -as.numeric(1/(1+t(phi) %*% phi)^(3/2))*phi
  a2 = as.numeric(1/sqrt(1+t(phi) %*% phi))*diag(1, 2) - as.numeric(1/(1+t(phi) %*% phi)^(3/2))*(phi %*% t(phi))
  return(matrix(rbind(a1, a2), 3, 2, byrow = T))
}

# Shortest arc-length distance
cdist = function(a, b){
  delta = abs(a - b)     
  pmin(delta, 2*pi - delta)
}


# Example data (angles in radians)
d = read.csv("dat.csv")  # Reading "dat.csv" after saving it in the current working directory

x0 = cbind(d[, 2], sin(d[, 3]), cos(d[, 3])) # Covariates
x = apply(x0, 2, function(x) x - mean(x))    # Centred covariates
y = d[, 1]                                   # Response



######### Fit Fisher-Lee model ###########

fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type = 'c-l'))



##########  Fit ZIvM 1 model ###########



int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

E.zivm = function(x, y, phi, gamma, a, mu, df = d.f){
  b = int_b(x_pr(beta(phi), x), df = df)
  delta = y - mu - 2*atan(as.vector(b %*% gamma))
  return(-sum((1 - a)*cos(delta)))
}


loglik.zivm1 = function(theta, x, y, df = d.f){
  q = ncol(x) - 1
  b = int_b(x_pr(beta(theta[4:(3+q)]), x), df = df)
  mu = theta[2] + 2*atan(as.vector(b %*% theta[-c(1:(3+q))]))
  sum(log(theta[1]*dvm(y, 0, 1000) + (1 - theta[1])*dvm(y, mu, theta[3])))
}


lc.zivm1 = function(x, y, init, df = d.f, tol = 1e-3, max.iter = 200){
  theta.0 = init
  diff = 1
  iter = 0
  q = ncol(x) - 1
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[4:(3+q)]), x), df = df)
    mu = theta.0[2] + 2*atan(as.vector(b %*% theta.0[-c(1:(3+q))]))
    a = (theta.0[1]*dvm(y, 0, 1000))/(theta.0[1]*dvm(y, 0, 1000) + (1 - theta.0[1])*dvm(y, mu, theta.0[3]))
    
    gamma_new = optim(theta.0[-c(1:(3+q))], function(gamma) E.zivm(x, y, theta.0[4:(3+q)], gamma, a, theta.0[2], df = df), method = "BFGS", control = list(maxit = 2000))$par
    phi_new = optim(theta.0[4:(3+q)], function(phi) E.zivm(x, y, phi, gamma_new, a, theta.0[2], df = df), method = "BFGS", control = list(maxit = 2000))$par
    
    p_new = mean(a)
    b = int_b(x_pr(beta(phi_new), x), df = df)
    del = y - 2*atan(as.vector(b %*% gamma_new))
    C = weighted.mean(cos(del), 1 - a)
    D = weighted.mean(sin(del), 1 - a)
    mu_new = atan2(D, C) 
    kappa_new = A1inv(sqrt(C^2 + D^2))
    theta.0 = c(p_new, mu_new, kappa_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)  
  } else {
    return(NULL)      
  }
}


zivm1 = function(x, y, n.knots = 3, m = 200){
  #Fit Fisher-Lee model
  fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type = 'c-l'))
  
  d.f = n.knots + 2
  # Initial values
  p0 = sum(y == 0)/length(y)
  beta0 = fl$coefficients
  phi0 = beta0[-1]/beta0[1]
  phi.init = matrix(phi0, nrow = m, ncol = length(phi0), byrow = TRUE)
  gamma.init = matrix(rnorm(m*d.f, 0, 0.3), nrow = m, ncol = d.f, byrow = TRUE)
  init.em = cbind(p0, as.numeric(fl$mu), fl$kappa, phi.init, gamma.init)
  
  # Parallel computation 
  ncores = detectCores()
  cl = makeCluster(ncores - 1)
  clusterEvalQ(cl, {library(splines2)
    library(CircStats)
    library(circular)
  })
  clusterExport(cl, c("l2", "beta", "x_pr", "int_b", "E.zivm", "lc.zivm1", "loglik.zivm1", "x", "y", "init.em", "d.f"),
                envir = environment())
  res = parLapply(cl, 1:m, function(i) {
    tryCatch({
      theta.em = lc.zivm1(x, y, init.em[i,], df = d.f)
      if (is.null(theta.em)) stop("EM did not converge")
      ll = loglik.zivm1(theta.em, x, y, df = d.f)
      c(theta.em, ll)
    }, error = function(e) NULL)
  })
  res = Filter(Negate(is.null), res)
  res = if (length(res) > 0) do.call(rbind, res) else NULL
  stopCluster(cl)
  
  ll = res[, ncol(res)]
  theta.em = res[, -ncol(res)]
  theta_hat = theta.em[which.max(ll),]
  aic = 2*length(theta_hat) - 2*loglik.zivm1(theta_hat, x, y, df = d.f)
  n = length(y)
  bic = length(theta_hat)*log(n) - 2*loglik.zivm1(theta_hat, x, y, df = d.f)
  q = ncol(x) - 1
  est = c(theta_hat[1:3], beta(theta_hat[4:(3+q)]), theta_hat[-c(1:(3+q))])
  cov = solve(hessian(function(theta) -loglik.zivm1(theta, x, y, df = d.f), theta_hat))
  se = sqrt(diag(cov))
  cov.beta = jac.beta(theta_hat[4:(3+q)]) %*% cov[4:(3+q), 4:(3+q)] %*% t(jac.beta(theta_hat[4:(3+q)]))
  se.beta = sqrt(diag(cov.beta))
  se_est = c(se[1:3], se.beta, se[-c(1:(3+q))])
  b = est[4:(4+q)]
  g = est[-c(1:(4+q))]
  
  # Estimated mean direction
  B_0 = int_b(x_pr(b, x), df = d.f)
  f_0 = as.vector(B_0 %*% g)
  mu_ZIVM1 = (est[2] + 2*atan(f_0)) %% (2*pi)
  
  # Sum of residuals (to be computed for non-zero responses and corresponding estimated mean directions)
  id = (y != 0)
  e1 = sum(cdist(y[id], mu_ZIVM1[id]))
  
  #CRPS
  n.sim = 1e6
  set.seed(1)
  sim = runif(n.sim, 0, 2*pi)
  sim1 = runif(n.sim, 0, 2*pi)
  sim2 = runif(n.sim, 0, 2*pi)
  crp11 = sapply(1:n, function(i){
    2*pi * mean(cdist(sim, y[i]) * (est[1]*(CircStats::dvm(sim, 0, 1000)) + (1 - est[1])*(CircStats::dvm(sim, mu_ZIVM1[i], est[3]))))
  })
  crp12 = sapply(1:n, function(i){
    4*pi^2 * mean(cdist(sim1, sim2) *
                    (est[1]*(CircStats::dvm(sim1, 0, 1000)) + (1 - est[1])*(CircStats::dvm(sim1, mu_ZIVM1[i], est[3]))) *
                    (est[1]*(CircStats::dvm(sim2, 0, 1000)) + (1 - est[1])*(CircStats::dvm(sim2, mu_ZIVM1[i], est[3]))))
  })
  crps1 = crp11 - 0.5*crp12
  
  result = list(est[1:(4+q)], se_est[1:(4+q)], aic, bic, e1, mean(crps1))
  names(result)=c("Estimate", "Std.Error", "AIC", "BIC", "Sum_Residual", "Avg.CRPS")
  names(result[[1]]) = names(result[[2]]) = c("p", "mu", "kappa", paste0("beta_", 1:(q+1)))
  result = lapply(result, round, 4)
  result$theta = theta_hat
  result$cov = cov
  result$d.f = d.f
  return(result)
}

t1 = Sys.time()
fit = zivm1(x, y)
t2 = Sys.time()
t2 - t1
fit



####### LOOCV ########



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
cv_zivm1 = numeric(n)
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
  
  fit = zivm1(x.train, y.train)
  theta = fit$theta
  b = beta(theta[4:(3+q)])
  g = theta[-c(1:(3+q))]
  
  xp.train = x_pr(b, x.train)
  xp.test  = x_pr(b, x.test)
  
  xp.test = pmin(pmax(xp.test, min(xp.train)), max(xp.train))
  
  if(length(unique(round(xp.train, 8))) <= 3){
    pred = NA
  } else {
    B_0 = try(int1_b(xp.test, xp.train, df = fit$d.f))
    if(inherits(B_0, "try-error")){
      pred = NA
    } else {
      f_0 = sum(g*B_0[1, ])
      pred = (theta[2] + 2*atan(f_0)) %% (2*pi)
    }
}
  if(all(is.na(pred))){
    cv_zivm1[i] = NA
  } else {
    crp1 = 2*pi * mean(cdist(sim, y.test)*(theta[1]*CircStats::dvm(sim, 0, 1000) +
                                              (1 - theta[1])*CircStats::dvm(sim, pred, theta[3])))
    
    crp2 = 4*pi^2 * mean(cdist(sim1, sim2) *(theta[1]*(CircStats::dvm(sim1, 0, 1000)) + (1 - theta[1])*(CircStats::dvm(sim1, pred, theta[3]))) *
                           (theta[1]*(CircStats::dvm(sim2, 0, 1000)) + (1 - theta[1])*(CircStats::dvm(sim2, pred, theta[3]))))
    
    crp = crp1 - 0.5*crp2
    cv_zivm1[i] = mean(crp)
   }
 }
})
t2 = Sys.time()
t2 - t1


# Average CV error
round(mean(cv_zivm1), 4)



############################################



# Plotting estimated function

theta_hat = fit$theta
cov = fit$cov

q = ncol(x) - 1
b = beta(theta_hat[4:(3+q)])
g = theta_hat[-c(1:(3+q))]

N = 100
u1 = seq(min(x_pr(b, x)), max(x_pr(b, x)), length = N)
B = int1_b(u1, x_pr(b, x), df = fit$d.f)

V = cov[-c(1:(3+q)), -c(1:(3+q))]
fhat = as.vector(B %*% g)
var.h = rowSums((B %*% V) * B)
lower = fhat - 1.96 * sqrt(var.h)
upper = fhat + 1.96 * sqrt(var.h)

#N=3
par(mar = c(5, 4.5, 4, 2) + 0.1)
plot(u1, fhat, type = "l", lwd = 2, xlab = expression(hat(beta)^T * x), ylab = expression(hat(h)(hat(beta)^T * x)),
     ylim = range(c(fhat, lower, upper), na.rm = TRUE))
lines(u1, lower, lty = 2)
lines(u1, upper, lty = 2)
h = function(x) x
curve(h, add = TRUE, lty = 4, lwd = 2)
legend("bottomleft", legend = c("ZIvM 1", "Confidence band for ZIvM 1", "Fisher-Lee"), lty = c(1, 2, 6), lwd = c(2, 1, 2), cex = 1.1, y.intersp = 0.7)   


# Plotting Donut-plots of ZIvM 1 (N=3) and Fisher-Lee

mu_fl = (as.numeric(fl$mu) + 2*atan(as.vector(x %*% fl$coefficients))) %% (2*pi) # Estimated mean direction in Fisher-Lee model

n.knots = 3
d.f = n.knots + 2
B_0 = int_b(x_pr(b, x), df = d.f)
f_0 = as.vector(B_0 %*% g)
mu_ZIVM1 = (theta_hat[2] + 2*atan(f_0)) %% (2*pi)


plot_circ_diag = function(y, y_hat) {
  stopifnot(length(y) == length(y_hat))
  delta = (y_hat - y) %% (2*pi) 
  r = 1 + cos(delta)
  x = r * cos(y_hat)
  y_cart = r * sin(y_hat)
  plot(0, 0, type = "n", asp = 1, xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5), axes = FALSE, xlab = "", ylab = "")
  symbols(rep(0, 2), rep(0, 2), circles = c(1, 2), inches = FALSE, add = TRUE, lwd = 1)
  points(0, 0, pch = 19, col = "black")
  
  cw  = delta > pi                # Clockwise 
  acw = delta > 0 & delta <= pi   # Anticlockwise 
  
  points(x[cw],  y_cart[cw],  pch = 19, col = "blue", cex = 2)   # Filled
  points(x[acw], y_cart[acw], pch = 1,  col = "blue", cex = 2)   # Open
}

par(mfrow = c(1, 2), mar = c(3, 2, 2, 1), oma = c(0, 0, 0, 0))    # shrink subplot margins
id = (y != 0)
plot_circ_diag(y[id], mu_fl[id])
plot_circ_diag(y[id], mu_ZIVM1[id])














