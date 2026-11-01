library(circular)
library(CircStats)
library(Directional)
library(parallel)


l2 = function(v){
  sqrt(sum(v^2))
}

cdist = function(a, b){
  delta = abs(a - b)     
  pmin(delta, 2*pi - delta)
}

# SPML density
fspml = function(theta, gamma, omega){
  z = gamma*cos(theta - omega)
  (1/(2*pi))*exp(-gamma^2/2) *
    (1 + z * pnorm(z) / dnorm(z))
}


# Example data (angles in radians)
d = read.csv("dat.csv")  # Reading "dat.csv" after saving it in the current working directory

x0 = cbind(d[, 2], sin(d[, 3]), cos(d[, 3])) # Covariates
x = apply(x0, 2, function(x) x - mean(x))    # Centred covariates
y = d[, 1]                                   # Response

n = length(y)

id = (y != 0)  # Non-zero responses



###############  Circular CRPS and Sum of Residuals  ##################


circ.crps = function(y, dens, n.sim = 1e6){
  set.seed(1)
  sim = runif(n.sim, 0, 2*pi)
  sim1 = runif(n.sim, 0, 2*pi)
  sim2 = runif(n.sim, 0, 2*pi)
  
  crp1 = sapply(seq_along(y), function(i)
    2*pi*mean(cdist(sim, y[i])*dens(sim, i)))
  
  crp2 = sapply(seq_along(y), function(i)
    4*pi^2*mean(cdist(sim1, sim2)*dens(sim1, i)*dens(sim2, i)))
  
  round(mean(crp1 - 0.5*crp2), 4)
}

circ.res = function(y, mu){round(sum(cdist(y, mu)), 4)}


# Fisher-Lee

fl = suppressWarnings(lm.circular(y, x, init = rnorm(ncol(x), 0, 0.3), type = 'c-l')) 
mu_fl = (as.numeric(fl$mu) + 2*atan(as.vector(x %*% fl$coefficients))) %% (2*pi)
crps.fl = circ.crps(y, function(t, i) CircStats::dvm(t, mu_fl[i], fl$kappa))

E.fl = circ.res(y[id], mu_fl[id])  # Sum of residuals


# SPML

bML = round(spml.reg(y, x0, rads = TRUE, xnew = NULL, tol = 1e-06)$be, 4)
mu1 = cbind(1, x0) %*% bML
gamma = apply(mu1, 1, l2)
omega = atan2(mu1[,2], mu1[,1]) %% (2*pi) 
crps.spml = circ.crps(y, function(t, i) fspml(t, gamma[i], omega[i]))

E.spml = circ.res(y[id], omega[id])  # Sum of residuals


# CIPC

bML2 = round(cipc.reg(y, x0, rads = TRUE, xnew = NULL, tol = 1e-06)$be, 4)
mu2 = cbind(1, x0) %*% bML2
gamma1 = apply(mu2, 1, l2)
omega1 = atan2(mu2[,2], mu2[,1]) %% (2*pi)
crps.cipc = circ.crps(y, function(t, i) Directional::dcipc(t, omega1[i], gamma1[i], rads = TRUE))

E.cipc = circ.res(y[id], omega1[id])  # Sum of residuals

# GCPC

fit.gcpc = gcpc.reg(y, x0, rads = TRUE, xnew = NULL)
bML3 = round(fit.gcpc$be,4)
rho = fit.gcpc$rho
mu3 = cbind(1, x0) %*% bML3
gamma2 = apply(mu3, 1, l2)
omega2 = atan2(mu3[, 2], mu3[, 1]) %% (2*pi)
crps.gcpc = circ.crps(y, function(t, i) Directional::dgcpc(t, omega2[i], gamma2[i], rho, rads = TRUE))

E.gcpc = circ.res(y[id], omega2[id])  # Sum of residuals



# Average circular CRPS
c(FL = crps.fl, SPML = crps.spml, CIPC = crps.cipc, GCPC = crps.gcpc)

# Sum of residuals
c(FL = E.fl, SPML = E.spml, CIPC = E.cipc, GCPC = E.gcpc)




###############   LOOCV  ##################



n.sim = 1e5
set.seed(1)
sim = runif(n.sim, 0, 2*pi)
sim1 = runif(n.sim, 0, 2*pi)
sim2 = runif(n.sim, 0, 2*pi)

d12 = cdist(sim1, sim2)

cv.fun = function(i){
  
  train.idx = (seq_len(n) != i)
  
  # Split predictors and response
  x.train = x[train.idx, , drop = FALSE]
  x.test  = x[i, , drop = FALSE]
  
  x0.train = x0[train.idx, , drop = FALSE ]
  x0.test = x0[i, , drop = FALSE ]
  
  y.train = y[train.idx]
  y.test  = y[i]
  
  fl = suppressWarnings(lm.circular(y.train, x.train, init = rnorm(ncol(x), 0, 0.3), type='c-l'))   # FL
  bML = spml.reg(y.train, x0.train, rads = TRUE)$be    # SPML
  bML2 = cipc.reg(y.train, x0.train, rads = TRUE)$be   # CIPC

  fit.gcpc = gcpc.reg(y.train, x0.train, rads = TRUE)  # GCPC
  bML3 = fit.gcpc$be  
  rho = fit.gcpc$rho
  
  
  # FL
  pred =  (as.numeric(fl$mu)+2*atan(as.vector(x.test %*% fl$coefficients))) %% (2*pi)
  
  crp11 = 2*pi * mean(cdist(sim, y.test) * (CircStats::dvm(sim, pred, fl$kappa)))
  crp12 = 4*pi^2 * mean(d12 * (CircStats::dvm(sim1, pred, fl$kappa)) * (CircStats::dvm(sim2, pred, fl$kappa)))
  crp1 = crp11-0.5*crp12
  
  # SPML
  mu1 = cbind(1, x0.test) %*% bML
  gamma = apply(mu1, 1, l2)
  pred1 = atan2(mu1[, 2], mu1[, 1]) %% (2*pi)
  
  crp21 = 2*pi * mean(cdist(sim, y.test) * fspml(sim, gamma, pred1))
  crp22 = 4*pi^2 * mean(d12 * fspml(sim1, gamma, pred1) * fspml(sim2, gamma, pred1))
  crp2 = crp21 - 0.5*crp22
  
  # CIPC
  mu2 = cbind(1, x0.test) %*% bML2
  gamma1 = apply(mu2, 1, l2)
  pred2 = atan2(mu2[,2], mu2[,1]) %% (2*pi)
  
  crp31 = 2*pi*mean(cdist(sim, y.test) * Directional::dcipc(sim, pred2, gamma1, rads = TRUE))
  crp32 = 4*pi^2*mean(d12 * Directional::dcipc(sim1, pred2, gamma1, rads = TRUE) * Directional::dcipc(sim2, pred2, gamma1, rads = TRUE))
  crp3 = crp31 - 0.5*crp32
  
  # GCPC
  mu3 = cbind(1, x0.test) %*% bML3
  gamma2 = apply(mu3, 1, l2)
  pred3 = atan2(mu3[,2], mu3[,1]) %% (2*pi)
  
  crp41 = 2*pi*mean(cdist(sim, y.test) * Directional::dgcpc(sim, pred3, gamma2, rho, rads = TRUE))
  crp42 = 4*pi^2*mean(d12 * Directional::dgcpc(sim1, pred3, gamma2, rho, rads = TRUE) * Directional::dgcpc(sim2, pred3, gamma2, rho, rads = TRUE))
  crp4 = crp41 - 0.5*crp42
  
  c(crp1, crp2, crp3, crp4)
}

# Parallelisation
ncores = detectCores()
cl = makeCluster(ncores - 1)

clusterExport(cl, c("l2", "x", "x0", "y", "n", "n.sim",
                    "sim", "sim1", "sim2", "d12", "cdist", "fspml", "cv.fun"), envir = environment())

invisible(clusterEvalQ(cl, {
  library(circular)
  library(CircStats)
  library(Directional)
}))

t1 = Sys.time()
res = parSapply(cl, 1:n, cv.fun)
stopCluster(cl)
t2 = Sys.time()
t2 - t1

cv.fl = res[1, ]
cv.spml = res[2, ]
cv.cipc = res[3, ]
cv.gcpc = res[4, ]


# Average CV error
cv.error = cbind(cv.fl, cv.spml, cv.cipc, cv.gcpc)
round(apply(cv.error, 2, mean), 4)














