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

# Beta reparameterization
beta = function(phi) c(1/sqrt(1+t(phi)%*%phi), as.numeric(1/sqrt(1+t(phi) %*% phi))*phi)

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

p = function(delta, x){
  eta = as.vector(x %*% delta)
  1/(1+exp(-eta))
}

w = function(delta, x){
  diag(p(delta, x)*(1 - p(delta, x)))
}


########## ZIvM 2 ###########



# Objective function for updating (phi,gamma) in ZIvM models
E.zivm = function(x, y, phi, gamma, a, mu){
  b = int_b(x_pr(beta(phi), x))
  delta = y - mu - 2*atan(as.vector(b %*% gamma))
  return(-sum((1 - a)*cos(delta)))
}

# Log-likelihood for ZIvM 2 model
loglik.zivm2 = function(theta, x, z, y){
  q1 = ncol(x) - 1
  q2 = ncol(z)
  b = int_b(x_pr(beta(theta[(q2+3):(q2+2+q1)]), x))
  mu = theta[(q2+1)] + 2*atan(as.vector(b %*% theta[-c(1:(q2+2+q1))]))
  sum(log(p(theta[1:q2], z)*dvm(y, 0, 1000) + (1 - p(theta[1:q2], z))*dvm(y, mu, theta[(q2+2)]) ))
}

# Function returning estimates under the ZIvM 2 model
lc.zivm2 = function(x, z, y, init, tol = 1e-3, max.iter = 200){
  theta.0 = init  # initial for delta,mu,kappa,phi,gamma
  diff = 1
  iter = 0
  q1 = ncol(x) - 1
  q2 = ncol(z)
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[(q2+3):(q2+2+q1)]), x))
    mu = theta.0[(q2+1)] + 2*atan(as.vector(b %*% theta.0[-c(1:(q2+2+q1))]))
    alpha = p(theta.0[1:q2], z)
    a = (alpha*dvm(y, 0, 1000))/(alpha*dvm(y, 0, 1000)+(1 - alpha)*dvm(y, mu, theta.0[(q2+2)]))
    
    init1 = theta.0[(q2+3):(q2+2+q1)]  # initial for phi
    init2 = theta.0[-c(1:(q2+2+q1))]  # initial for gamma
    gamma_new = optim(init2, function(gamma) E.zivm(x, y, init1, gamma, a, theta.0[(q2+1)]), method = "BFGS", control = list(maxit = 2000))$par
    phi_new = optim(init1, function(phi) E.zivm(x, y, phi, gamma_new, a, theta.0[(q2+1)]), method = "BFGS", control = list(maxit = 2000))$par
    
    delta_new = theta.0[1:q2] + solve(t(z) %*% w(theta.0[1:q2], z) %*% z) %*% t(z) %*% (a - p(theta.0[1:q2], z))
    delta_new = as.numeric(delta_new)
    
    b = int_b(x_pr(beta(phi_new), x))
    del = y - 2*atan(as.vector(b %*% gamma_new))
    C = weighted.mean(cos(del), 1 - a)
    D = weighted.mean(sin(del), 1 - a)
    mu_new = atan2(D, C) 
    kappa_new = A1inv(sqrt(C^2 + D^2))
    theta.0 = c(delta_new, mu_new, kappa_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)      
  }
}

# Function returning estimates under the ZIvM 2 model given multiple random initial values                    
zivm2.em = function(x, z, y, init){
  M = nrow(init)
  theta.em = t(sapply(1:M, function(i) lc.zivm2(x, z, y, init[i,])))
  ll = sapply(1:M, function(i) loglik.zivm2(theta.em[i,], x, z, y))
  theta_hat = theta.em[which.max(ll),]
  return(theta_hat)
}                


########## ZIJP 2 #############


# Normalizing constant of Jones-Pewsey distribution
JPNCon = function(kappa, psi) {
  if (kappa < 0.001) {ncon = 1/(2*pi) ; return(ncon) }
  else {
    eps = 1e-6
    if (abs(psi) <= eps) { ncon = 1/(2*pi*I.0(kappa)) ; return(ncon) }
    else {
      intgrnd = function(x) { (cosh(kappa*psi)+sinh(kappa*psi)*cos(x))^(1/psi) }
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
  if (abs(kappa*psi) > 10) return(Inf)
  else { ncon = JPNCon(kappa, psi)
  return(-sum((1 - a)*log(djp(y, mu, kappa, psi, ncon)))) }  
}

# Objective function for updating (phi,gamma) in ZIJP models
E.zijp = function(x, y, phi, gamma, a, mu, kappa, psi){
  if(abs(kappa*psi) > 10) return(Inf)
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

# Log-likelihood for ZIJP 2 model
loglik.zijp2 = function(theta, x, z, y){
  q1 = ncol(x) - 1
  q2 = ncol(z)
  b = int_b(x_pr(beta(theta[(q2+4):(q2+3+q1)]), x))
  mu = theta[(q2+1)] + 2*atan(as.vector(b %*% theta[-c(1:(q2+3+q1))]))
  if(abs(theta[(q2+2)]*theta[(q2+3)]) > 10) return(-Inf)
  else {
    ncon = JPNCon(theta[(q2+2)], theta[(q2+3)])
    sum(log(p(theta[1:q2], z)*dvm(y, 0, 1000) + (1 - p(theta[1:q2], z))*djp(y, mu, theta[(q2+2)], theta[(q2+3)], ncon)))
  }
}

# Function returning estimates under the ZIJP 2 model
lc.zijp2 = function(x, z, y, init, tol = 1e-3, max.iter = 200){
  theta.0 = init
  diff = 1
  iter = 0
  q1 = ncol(x) - 1
  q2 = ncol(z)
  while(diff > tol && iter < max.iter){
    theta.1 = theta.0
    b = int_b(x_pr(beta(theta.0[(q2+4):(q2+3+q1)]), x))
    mu = theta.0[(q2+1)] + 2*atan(as.vector(b %*% theta.0[-c(1:(q2+3+q1))]))
    ncon = JPNCon(theta.0[(q2+2)], theta.0[(q2+3)])
    alpha = p(theta.0[1:q2], z)
    a = (alpha*dvm(y, 0, 1000))/(alpha*dvm(y, 0, 1000)+(1-alpha)*djp(y, mu, theta.0[(q2+2)], theta.0[(q2+3)], ncon))
    
    init1 = theta.0[(q2+4):(q2+3+q1)]  # initial for phi
    init2 = theta.0[-c(1:(q2+3+q1))]   # initial for gamma
    gamma_new = optim(init2, function(gamma) E.zijp(x, y, init1, gamma, a, theta.0[(q2+1)], theta.0[(q2+2)], theta.0[(q2+3)]), method = "BFGS", control = list(maxit = 2000))$par
    phi_new = optim(init1, function(phi) E.zijp(x, y, phi, gamma_new, a, theta.0[(q2+1)], theta.0[(q2+2)], theta.0[(q2+3)]), method = "BFGS", control = list(maxit = 2000))$par
    
    delta_new = theta.0[1:q2] + solve(t(z) %*% w(theta.0[1:q2], z) %*% z) %*% t(z) %*% (a - p(theta.0[1:q2], z))
    delta_new = as.numeric(delta_new)
    
    out = optim(theta.0[(q2+1):(q2+3)], function(p) zijpnll(x, y, a, p, phi_new, gamma_new), method = "L-BFGS-B", lower = c(-pi, 0, -Inf), upper = c(pi, Inf, Inf), control = list(maxit = 2000))
    mu_new = out$par[1] ; kappa_new = out$par[2] ; psi_new = out$par[3]
    theta.0 = c(delta_new, mu_new, kappa_new, psi_new, phi_new, gamma_new)
    diff = l2(theta.0 - theta.1)
    iter = iter + 1
  }
  if (diff <= tol) {
    return(theta.0)   
  } else {
    return(NULL)      
  }
}

# Function returning estimates under the ZIJP 2 model given multiple random initial values
zijp2.em = function(x, z, y, init){
  M = nrow(init)
  theta.list = vector("list", M)
  ll = rep(-Inf, M)
  for(i in 1:M){
    fit = lc.zijp2(x, z, y, init[i,])
    if(!is.null(fit)){
      theta.list[[i]] = fit
      ll[i] = loglik.zijp2(fit, x, z, y)
    }
  }
  if(all(is.infinite(ll))){
    stop("No initial value converged")
  }
  theta_hat = theta.list[[which.max(ll)]]
  return(theta_hat)
}


####### Simulation #######


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
delta0 = -1; delta1 = -0.5; delta2 = 0.5
delta = c(delta0, delta1, delta2)
mu0 = pi/3
kappa = 1.5
psi = -0.75
ncon = JPNCon(kappa, psi)
beta.0 = c(1, 1)/sqrt(2)
phi = beta.0[-1]/beta.0[1]

f0 = function(x) 1/(1 + exp(-3*(x - 0.2)))
f_0 = function(x) f0(x) - f0(0)             # True single-index function

# Random initial values for spline coefficients
m = 20

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
  "l2", "x_pr", "beta", "p", "w", "int_b", "int1_b",
  
  # ZIvM 2 functions
  "E.zivm", "loglik.zivm2", "lc.zivm2", "zivm2.em",
  
  # ZIJP 2 functions
  "JPNCon", "djp", "jpsim", "zijpnll", "E.zijp", "loglik.zijp2", "lc.zijp2", "zijp2.em",
  
  # Simulation parameters and objects
  "delta0", "delta1", "delta2", "delta", "mu0", "kappa", "psi", "ncon", "beta.0", "phi", "f0", "f_0", "m", "d.f", "t"
), envir = environment())

res = parLapply(cl, 1:550, function(j) {
  tryCatch({
    # Data generation under ZIJP 2
    set.seed(j)
    n = 200
    x1 = rnorm(n, 0, 0.6)
    x2 = runif(n, -1.1, 1.1)
    x = cbind(x1, x2)
    z = cbind(1, x)
    U = runif(n)
    y = rep(0, n)
    mu = mu0 + 2*atan(as.vector(f_0(x %*% beta.0)))
    prob = p(delta, z)
    for (i in 1:n) {
      if (U[i] < prob[i]) {
        y[i] = rvm(1, 0, 1000)
      } else {
        y[i] = jpsim(1, mu[i], kappa, psi, ncon)
      }
    }
    gamma.init = matrix(rnorm(m*d.f, 0, 0.3), m, d.f, byrow = TRUE)
    init1 = cbind(delta0, delta1, delta2, mu0, kappa, psi, phi, gamma.init)      # Initial values for ZIJP 2
    fit = zijp2.em(x, z, y, init1)                                               # ZIJP 2 fit
    
    #init2 = cbind(delta0, delta1, delta2, mu0, kappa, phi, gamma.init)           # Initial values for ZIvM 2
    #fit = zivm2.em(x, z, y, init2)                                               # ZIvM 2 fit
    
    q1 = ncol(x) - 1
    q2 = ncol(z)
    
    ## Estimated beta
    b = unname(beta(fit[(q2+4):(q2+3+q1)]))          # ZIJP 2
    #b = unname(beta(fit[(q2+3):(q2+2+q1)]))         # ZIvM 2
    
    ## Projected training values
    eta = x_pr(b, x)
    
    ## Evaluate fitted basis on common grid
    B = int1_b(t, eta)
    
    ## Estimated spline coefficients
    g = fit[-c(1:(q2+3+q1))]         # ZIJP 2
    #g = fit[-c(1:(q2+2+q1))]        # ZIvM 2
    
    ## Fitted function
    hhat = as.vector(B %*% g)
    
    mu.true = (mu0 + 2*atan(f_0(t))) %% (2*pi)
    mu.hat  = (fit[q2+1] + 2*atan(hhat)) %% (2*pi)
    cmspe = mean(1 - cos(mu.hat - mu.true))            # Circular MSPE
    
    c(fit, cmspe, hhat)
    
  }, error = function(e) NULL)
})
stopCluster(cl)
t2 = Sys.time()
t2 - t1

res = Filter(Negate(is.null), res)
res = if (length(res) > 0) do.call(rbind, res) else NULL

res1 = res[seq_len(min(500, nrow(res))), ]    # Considering first 500 successful replications
q1 = 1   # ncol(x) - 1
q2 = 3   # ncol(z)

# Average circular MSPE 
p1 = q2+3+q1          # ZIJP 2
#p1 = q2+2+q1         # ZIvM 2
cMSPE = mean(res1[, p1+d.f+1])


### Plotting the average estimated single-index function with 95% pointwise confidence interval ###
hhat = res1[, -c(1:(p1+d.f+1))]
mean.fun = colMeans(hhat)
par(mar = c(5, 4.5, 4, 2) + 0.1)
lower = apply(hhat, 2, quantile, probs = 0.025)
upper = apply(hhat, 2, quantile, probs = 0.975)
plot(t, mean.fun, type = "l", xlab = expression(hat(beta)^T*x), ylab = expression(bar(h)(hat(beta)^T*x)), lty = 4, lwd = 2,
     ylim = range(c(mean.fun, lower, upper), na.rm = TRUE))
curve(f_0, add = T, lwd = 2)
lines(t, lower, lty = 2)
lines(t, upper, lty = 2)


# Parameter estimates (with their standard deviations) and average circular MSPE
#ZIJP 2
beta_hat = unname(t(sapply(res1[, (q2+4):(q2+3+q1)], beta))) 
res1 = cbind(res1[, 1:(q2+3)], beta_hat, res1[, -c(1:(q2+3+q1))])
sim = list(Estimate = apply(res1[, 1:(q2+4+q1)], 2, mean),
           Std.Error = apply(res1[, 1:(q2+4+q1)], 2, sd), cMSPE = cMSPE)

#ZIvM 2
beta_hat = unname(t(sapply(res1[, (q2+3):(q2+2+q1)], beta)))         
res1 = cbind(res1[, 1:(q2+2)], beta_hat, res1[, -c(1:(q2+2+q1))]) 
sim = list(Estimate = apply(res1[, 1:(q2+3+q1)], 2, mean),
           Std.Error = apply(res1[, 1:(q2+3+q1)], 2, sd), cMSPE = cMSPE)

sim$Estimate[q2+1] = circ.mean(res1[, q2+1])
sim$Std.Error[q2+1] = sd.circular(res1[, q2+1])
names(sim$Estimate) = names(sim$Std.Error) = c(paste0("delta_", 0:(q2-1)), "mu", "kappa", "psi", paste0("beta_", 1:(q1+1)))     # ZIJP 2
#names(sim$Estimate) = names(sim$Std.Error) = c(paste0("delta_", 0:(q2-1)), "mu", "kappa", paste0("beta_", 1:(q1+1)))            # ZIvM 2
lapply(sim, round, 4)



















































