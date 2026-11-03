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

# Beta reparameterisation
beta = function(phi) c(1/sqrt(1+t(phi) %*% phi), as.numeric(1/sqrt(1+t(phi) %*% phi))*phi)

# Jacobian matrix of beta with respect to phi
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


########## ZIvM 1 ###########

# Objective function for updating (phi,gamma) in ZIvM models
E.zivm = function(x, y, phi, gamma, a, mu){
  b = int_b(x_pr(beta(phi), x))
  delta = y - mu - 2*atan(as.vector(b %*% gamma))
  return(-sum((1 - a)*cos(delta)))
}

# Log-likelihood for ZIvM 1 model
loglik.zivm1 = function(theta, x, y){
  q = ncol(x) - 1
  b = int_b(x_pr(beta(theta[4:(3+q)]), x))
  mu = theta[2] + 2*atan(as.vector(b %*% theta[-c(1:(3+q))]))
  sum(log(theta[1]*dvm(y, 0, 1000) + (1 - theta[1])*dvm(y, mu, theta[3])))
}

# Function returning estimates under the ZIvM 1 model
lc.zivm1 = function(x, y, init, tol = 1e-3, max.iter = 200){
  theta.0 = init
  diff = 1
  iter = 0
  q = ncol(x) - 1
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[4:(3+q)]), x))
    mu = theta.0[2] + 2*atan(as.vector(b %*% theta.0[-c(1:(3+q))]))
    a = (theta.0[1]*dvm(y, 0, 1000))/(theta.0[1]*dvm(y, 0, 1000)+(1 - theta.0[1])*dvm(y, mu, theta.0[3]))
    
    gamma_new = optim(theta.0[-c(1:(3+q))], function(gamma) E.zivm(x, y, theta.0[4:(3+q)], gamma, a, theta.0[2]), method = "BFGS", control = list(maxit = 2000))$par
    phi_new = optim(theta.0[4:(3+q)], function(phi) E.zivm(x, y, phi, gamma_new, a, theta.0[2]), method = "BFGS", control = list(maxit = 2000))$par
    
    p_new = mean(a)
    b = int_b(x_pr(beta(phi_new), x))
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

# Function returning estimates under the ZIvM 1 model given multiple random initial values                    
zivm1.em = function(x, y, init){
  M = nrow(init)
  theta.em = t(sapply(1:M, function(i) lc.zivm1(x, y, init[i,])))
  ll = sapply(1:M, function(i) loglik.zivm1(theta.em[i,], x, y))
  theta_hat = theta.em[which.max(ll),]
  return(theta_hat)
}                     


########## ZIJP 1 #############


# Normalizing constant of Jones-Pewsey distribution
JPNCon = function(kappa, psi) {
  if (kappa < 0.001) {ncon = 1/(2*pi) ; return(ncon) }
  else {
    eps = 1e-6
    if (abs(psi) <= eps) { ncon = 1/(2*pi*I.0(kappa)) ; return(ncon) }
    else {
      intgrnd = function(x) { (cosh(kappa*psi) + sinh(kappa*psi)*cos(x))^(1/psi) }
      ncon = 1/integrate(intgrnd, lower = -pi, upper = pi)$value
      return(ncon) } }
}

# Density of Jones-Pewsey distribution
djp = function(theta, mu, kappa, psi, ncon) {
  if (kappa < 0.001) {pdfval = 1/(2*pi) ; return(pdfval)}
  else {
    eps = 1e-6
    if (abs(psi) <= eps) {
      pdfval = ncon*exp(kappa*cos(theta - mu)) ; return(pdfval) }
    else {
      pdfval = (cosh(kappa*psi)+sinh(kappa*psi)*cos(theta - mu))^(1/psi)
      pdfval = ncon*pdfval ; return(pdfval) } }
}

# Random number generation from Jones-Pewsey distribution
jpsim = function(n, mu, kappa, psi, ncon) {
  fmax = djp(mu, mu, kappa, psi, ncon) ; theta = 0
  for (j in 1:n) { stopgo = 0
  while (stopgo == 0) {
    u1 = runif(1, 0, 2*pi) ; pdfu1 = djp(u1, mu, kappa, psi, ncon)
    u2 = runif(1, 0, fmax)
    if (u2 <= pdfu1) { theta[j] = u1 ; stopgo = 1 }
  } }
  return(theta)
}

# Negative log-likelihood of Zero-inflated Jones-Pewsey (ZIJP) distribution
zijpnll = function(x, y, a, p, phi, gamma) {
  b = int_b(x_pr(beta(phi), x))
  mu = p[1] + 2*atan(as.vector(b %*% gamma))
  kappa = p[2] ; psi = p[3]
  if (abs(kappa*psi) > 10) return(9999)
  else { ncon = JPNCon(kappa, psi)
  return(-sum((1 - a)*log(djp(y, mu, kappa, psi, ncon)))) }  
}

# Objective function for updating (phi,gamma) in ZIJP models
E.zijp = function(x, y, phi, gamma, a, mu, kappa, psi){
  if(abs(kappa*psi) > 10) return(9999)
  b = int_b(x_pr(beta(phi), x))
  delta = y - mu - 2*atan(as.vector(b %*% gamma))
  eps = 1e-6
  if(abs(psi) <= eps){
    return(-sum((1 - a)*cos(delta)))
  } else {
    S = sum((1 - a)*log(cosh(kappa*psi) + sinh(kappa*psi)*cos(delta)))
    return(-(1/psi)*S)
  }
}

# Log-likelihood for ZIJP 1 model
loglik.zijp1 = function(theta, x, y){
  q = ncol(x) - 1
  b = int_b(x_pr(beta(theta[5:(4+q)]), x))
  mu = theta[2] + 2*atan(as.vector(b %*% theta[-c(1:(4+q))]))
  if(abs(theta[3]*theta[4]) > 10) return(-9999)
  else {
    ncon = JPNCon(theta[3], theta[4])
    sum(log(theta[1]*dvm(y, 0, 1000) + (1 - theta[1])*djp(y, mu, theta[3], theta[4], ncon)))
  }
}

# Function returning estimates under the ZIJP 1 model
lc.zijp1 = function(x, y, init, tol = 1e-3, max.iter = 200){
  theta.0 = init
  diff = 1
  iter = 0
  q = ncol(x) - 1
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[5:(4+q)]), x))
    mu = theta.0[2] + 2*atan(as.vector(b %*% theta.0[-c(1:(4+q))]))
    ncon = JPNCon(theta.0[3], theta.0[4])
    a = (theta.0[1]*dvm(y, 0, 1000))/(theta.0[1]*dvm(y, 0, 1000) + (1 - theta.0[1])*djp(y, mu, theta.0[3], theta.0[4], ncon))
    
    gamma_new = optim(theta.0[-c(1:(4+q))], function(gamma) E.zijp(x, y, theta.0[5:(4+q)], gamma, a, theta.0[2], theta.0[3], theta.0[4]), method = "BFGS", control=list(maxit = 2000))$par
    phi_new = optim(theta.0[5:(4+q)], function(phi) E.zijp(x, y, phi, gamma_new, a, theta.0[2], theta.0[3], theta.0[4]), method = "BFGS", control = list(maxit = 2000))$par
    out = optim(theta.0[2:4], function(p) zijpnll(x, y, a, p, phi_new, gamma_new), method = "L-BFGS-B", lower = c(-pi, 0, -Inf), upper = c(pi, Inf, Inf), control = list(maxit = 2000))
    mu_new = out$par[1] ; kappa_new = out$par[2] ; psi_new = out$par[3]
    p_new = mean(a)
    theta.0 = c(p_new, mu_new, kappa_new, psi_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)      
  }
}

# Function returning estimates under the ZIJP 1 model given multiple random initial values
zijp1.em = function(x, y, init){
  M = nrow(init)
  theta.list = vector("list", M)
  ll = rep(-Inf, M)
  for(i in 1:M){
    fit = lc.zijp1(x, y, init[i,])
    if(!is.null(fit)){
      theta.list[[i]] = fit
      ll[i] = loglik.zijp1(fit, x, y)
    }
  }
  if(all(is.infinite(ll))){
    stop("No initial value converged")
  }
  theta_hat = theta.list[[which.max(ll)]]
  return(theta_hat)
}
                


############ Simulation #############



n.knots = 3               # Number of interior knots
d.f = n.knots + 2         # Degrees of freedom of the spline function

# Integrated B-spline basis function               
int_b = function(x, df = d.f, degree = 2) {
  ibs_x = ibs(x, df, degree = degree)
  ibs_0 = predict(ibs_x, 0)
  sweep(ibs_x, 2, ibs_0, "-")
}

int1_b = function(x.new, x.ref, df = d.f, degree = 2){
  ibs_ref = ibs(x.ref, df = df, degree = degree)
  ibs_new = predict(ibs_ref, x.new)
  ibs_0   = predict(ibs_ref, 0)
  sweep(ibs_new, 2, ibs_0, "-")
}

#Equispaced points between -1 and 1 for plotting the average estimated single-index function               
N = 100
t = seq(-1, 1, length = N)

# True values of parameters
p = 0.1
mu0 = pi/5
kappa = 2
psi = -0.5
ncon = JPNCon(kappa, psi)
beta.0 = c(1, 0)  
#beta.0 = c(1/2, sqrt(3)/2)
#beta.0 = c(1/sqrt(3), sqrt(2/3))
phi = beta.0[-1]/beta.0[1]

# True single-index function
f_0 = function(t) sin(5*t)/2  
#f_0 = function(t) 2*t^2
#f_0 = function(t) t
                
m = 20  # Number of random initial values for spline coefficients

# Parallel computation of estimates over multiple replications under ZIJP 1 
ncores = detectCores()
t1 = Sys.time()
cl = makeCluster(ncores-1)

invisible(clusterEvalQ(cl, {
  library(circular)
  library(CircStats)
  library(splines2)
}))

clusterExport(cl, c(
  # Common functions
  "l2", "beta", "x_pr", "int_b", "int1_b",
  
  # ZIvM 1 functions
  "E.zivm", "loglik.zivm1", "lc.zivm1", "zivm1.em",
  
  # ZIJP 1 functions
  "JPNCon", "djp", "jpsim", "E.zijp", "zijpnll", "loglik.zijp1", "lc.zijp1", "zijp1.em",
  
  # Simulation parameters and objects
  "f_0", "mu0", "beta.0", "p", "kappa", "psi", "ncon", "phi", "m", "d.f", "t"
), envir = environment())

res = parLapply(cl, 1:1050, function(j) {
  tryCatch({
    # Data generation under ZIJP 1
    set.seed(j)
    n = 200
    x1 = rnorm(n, 0, 0.6)
    x2 = runif(n, -1.1, 1.1)
    x = cbind(x1, x2)
    U = runif(n)
    y = rep(0, n)
    mu = mu0 + 2*atan(as.vector(f_0(x %*% beta.0)))
    for (i in 1:n) {
      if (U[i] < p) {
        y[i] = rvm(1, 0, 1000)
      } else {
        y[i] = jpsim(1, mu[i], kappa, psi, ncon)
      }
    }
    gamma.init = matrix(rnorm(m*d.f, 0, 0.3), m, d.f)
    init1 = cbind(p, mu0, kappa, psi, phi, gamma.init)         # Initial values for ZIJP 1
    fit = zijp1.em(x, y, init1)                                # ZIJP 1 fit
    
    #init2 = cbind(p, mu0, kappa, phi, gamma.init)              # Initial values for ZIvM 1
    #fit = zivm1.em(x, y, init2)                                # ZIvM 1 fit
    
    q = ncol(x) - 1
    
    ## Estimated beta
    b = unname(beta(fit[5:(4+q)]))    # ZIJP 1
    #b = unname(beta(fit[4:(3+q)]))    # ZIvM 1
    
    ## Projected training values
    eta = x_pr(b, x)
    
    ## Evaluate fitted basis on common grid
    B = int1_b(t, eta)
    
    ## Estimated spline coefficients
    g = fit[-c(1:(4+q))]             # ZIJP 1
    #g = fit[-c(1:(3+q))]             # ZIvM 1
    
    ## Fitted function
    hhat = as.vector(B %*% g)
    
    mu.true = (mu0 + 2*atan(f_0(t))) %% (2*pi)
    mu.hat = (fit[2] + 2*atan(hhat)) %% (2*pi)
    cmspe = mean(1 - cos(mu.hat - mu.true))     # Circular MSPE
    
    c(fit, cmspe, hhat)
    
  }, error = function(e) NULL)
})
stopCluster(cl)
t2 = Sys.time()
t2 - t1

res = Filter(Negate(is.null), res)
res = if (length(res) > 0) do.call(rbind, res) else NULL

res1 = res[seq_len(min(1000, nrow(res))), ]    # Considering first 1000 successful replications
q = 1  # ncol(x) - 1

# Average circular MSPE 
cMSPE = mean(res1[, 5+q+d.f])   # ZIJP 1
#cMSPE = mean(res1[, 4+q+d.f])   # ZIvM 1


### Plotting the average estimated single-index function with 95% pointwise confidence interval ###

hhat = res1[, -c(1:(q+5+d.f))]     # ZIJP 1
#hhat = res1[, -c(1:(q+4+d.f))]     # ZIvM 1

mean.fun = colMeans(hhat)
par(mar = c(5, 4.5, 4, 2) + 0.1)
lower = apply(hhat, 2, quantile, probs = 0.025)
upper = apply(hhat, 2, quantile, probs = 0.975)
plot(t, mean.fun, type = "l", xlab = expression(hat(beta)^T*x), ylab = expression(bar(h)(hat(beta)^T*x)), lty = 4, lwd = 2, 
     ylim = range(c(mean.fun, lower, upper), na.rm = TRUE))
curve(f_0, add = T, lwd = 2)
lines(t, lower, lty = 2)
lines(t, upper, lty = 2)


### Parameter estimates (with their standard deviations) and average circular MSPE ###

# ZIJP 1
beta_hat = unname(t(sapply(res1[, 5:(4+q)], beta)))
res1 = cbind(res1[, 1:4], beta_hat, res1[, -c(1:(4+q))])
sim = list(Estimate = apply(res1[, 1:(5+q)], 2, mean),
           Std.Error = apply(res1[, 1:(5+q)], 2, sd), cMSPE = cMSPE)

# ZIvM 1
#beta_hat = unname(t(sapply(res1[, 4:(3+q)], beta)))
#res1 = cbind(res1[, 1:3], beta_hat, res1[, -c(1:(3+q))])
#sim = list(Estimate = apply(res1[, 1:(4+q)], 2, mean),
#           Std.Error = apply(res1[, 1:(4+q)], 2, sd), cMSPE = cMSPE)

sim$Estimate[2] = circ.mean(res1[, 2])
sim$Std.Error[2] = sd.circular(res1[, 2])
names(sim$Estimate) = names(sim$Std.Error) = c("p", "mu", "kappa", "psi", paste0("beta_", 1:(q+1)))   # ZIJP 1
#names(sim$Estimate) = names(sim$Std.Error) = c("p", "mu", "kappa", paste0("beta_", 1:(q+1)))         # ZIvM 1
lapply(sim, round, 4)













































