library(geoR)
library(gwselect)
library(doMC)
registerCores(n=3)

seeds = as.vector(read.csv("seeds.csv", header=FALSE)[,1])
B = 100
N = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
settings = 9
tau = rep(0, settings)
rho = rep(c(0,0.5,0.9), settings/3)
#sigma.tau = rep(0.05, settings)
#sigma = rep(c(0.5,1), settings/2)
params = data.frame(tau, rho)

#Read command-line parameters
cluster=NA
process=1

#Simulation parameters are based on the value of process
setting = process %/% B + 1
parameters = params[setting,]
set.seed(seeds[process+1])

parameters = list(tau=0, rho=0, sigma.tau=0, sigma=0.5)

#Generate the covariates:
if (parameters[['tau']] > 0) {
    d1 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d2 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d3 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d4 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
    d5 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))$data
} else {
    d1 = rnorm(N**2, mean=0, sd=1)
    d2 = rnorm(N**2, mean=0, sd=1)
    d3 = rnorm(N**2, mean=0, sd=1)
    d4 = rnorm(N**2, mean=0, sd=1)
    d5 = rnorm(N**2, mean=0, sd=1)
}

loc.x = rep(coord, times=N)
loc.y = rep(coord, each=N)

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1, d2, d3, d4, d5)) %*% L
    
#
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)

if ((setting-1) %/% 2 == 0) {
    B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
} else if ((setting-1) %/% 2 == 1) {
    B1 = matrix(rep(coord, N), N, N)
} else if ((setting-1) %/% 2 == 2) {
    Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
    Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
    D = (Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2
    d = D[435,]
    B1 = matrix(max(d)-d, N, N)
}


#if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N**2, mean=0, sd=parameters[['sigma']])}
#if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(parameters[['sigma']]**2,parameters[['sigma.tau']]))$data}



nu = X1*B1
Y = rpois(length(nu), lambda=exp(nu))

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
fitloc = cbind(rep(seq(0,1, length.out=N), each=N), rep(seq(0,1, length.out=N), times=N))

vars = cbind(B1=as.vector(B1!=0))#, B2=as.vector(B2!=0), B3=as.vector(B3!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}


#Find the optimal bandwidth and use it to generate a model:
bw = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='poisson', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="BIC", gweight=spherical, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE, AICc=TRUE)
model = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='poisson', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=100, mode.select='BIC', bw=bw[['bw']], gweight=spherical, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE, AICc=TRUE)
