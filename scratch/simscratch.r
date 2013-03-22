library(geoR)
library(gwselect)
library(doMC)
registerCores(n=3)

seeds = as.vector(read.csv("seeds.csv", header=FALSE)[,1])
B = 100
N = 30
N.full = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(c(0.03, 0.1), each=9)
rho = rep(rep(c(0, 0.5, 0.8), each=3), times=2)
sigma.tau = rep(c(0, 0.03, 0.1), times=6)
b = 25
B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)

params = data.frame(tau, rho, sigma.tau)

#Read command-line parameters
#args = commandArgs(trailingOnly=TRUE)
#cluster = as.integer(args[1])
#process = as.integer(args[2])
cluster=NA
process=2

#Simulation parameters are based on the value of process
setting = process %/% B + 1
parameters = params[setting,]

#Get two (independent) Gaussian random fields:
set.seed(seeds[process+1])
d1 = grf(n=N.full**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d2 = grf(n=N.full**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d3 = grf(n=N.full**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d4 = grf(n=N.full**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))
d5 = grf(n=N.full**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['tau']]))

loc.x = d1$coords[,1]
loc.y = d1$coords[,2]

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(parameters[['rho']], 5, 5)
diag(S) = rep(1, 5)
L = chol(S)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1$data, d2$data, d3$data, d4$data, d5$data)) %*% L
    
#
X1 = matrix(D[,1], N.full, N.full)
X2 = matrix(D[,2], N.full, N.full)
X3 = matrix(D[,3], N.full, N.full)
X4 = matrix(D[,4], N.full, N.full)
X5 = matrix(D[,5], N.full, N.full)

#if (parameters[['function.type']] == 'step') {B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N.full), N.full, N.full)}
#if (parameters[['function.type']] == 'gradient') {B1 = matrix(rep(1-coord, N.full), N.full, N.full)}


if (parameters[['sigma.tau']] == 0) {epsilon = rnorm(N.full**2, mean=0, sd=1)}
if (parameters[['sigma.tau']] > 0) {epsilon = grf(n=N.full**2, grid='reg', cov.model='exponential', cov.pars=c(1,parameters[['sigma.tau']]))$data}

#
mu = X1*B1
Y = mu + epsilon

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
fitloc = cbind(rep(seq(0,1, length.out=N), each=N), rep(seq(0,1, length.out=N), times=N))

vars = as.vector(B1!=0)
oracle = list()
for (i in 1:N**2) { 
    if (vars[i]) {
        oracle[[i]] = c("X1")
    } else { 
        oracle[[i]] = character(0)
    }
}

#Find the optimal bandwidth and use it to generate a model:
#bw = gwlars.sel(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,0.5), gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)
#model = gwlars(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)


#bw.glmnet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,0.5), gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)
#model.glmnet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.glmnet, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)

#bw.enet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,0.5), gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)
#model.enet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.enet, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)

#bw.oracular = gwlars.sel(Y~X1+X2+X3+X4+X5-1, data=sim, oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,0.5), gweight=bisquare, tol=0.01, method='dist', parallel=FALSE, interact=TRUE)
#model.oracular = gwlars(Y~X1+X2+X3+X4+X5-1, data=sim, oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.oracular, gweight=bisquare, tol=0.01, method='dist', simulation=TRUE, parallel=FALSE, interact=TRUE)


n=900
coords = sim[,c('loc.x','loc.y')]
Xmat = matrix(rep(coords[,1], times=n), n, n)
Ymat = matrix(rep(coords[,2], times=n), n, n)
D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
d = as.vector(D[280,])
bw = 0.25
w = bisquare(d, bw)
indx = which(w>0)

W = diag(sqrt(w[indx]))
X = as.matrix(sim[indx,c('X1','X2','X3','X4','X5')])
Y = as.matrix(sim$Y[indx])
wm = sum(w[indx]*Y) / sum(w)

WX = W %*% cbind(1,X)
WI = W %*% matrix(1, nrow=dim(X)[1], ncol=1)
WX.2 = W %*% X
WY = W %*% Y
WY.c = WY - mean(WY)
WX.c = scale(WX, scale=FALSE)
cm = apply(X, 2, function(x) {sum(x*w[indx])})/sum(w)
wcm = apply(WX[,2:6], 2, function(x) {sum(x*sqrt(w[indx]))})/sum(sqrt(w))

lsm = lsfit(x=X, y=Y, wt=w[indx])
lsm.2 = lm(Y~X, weights=w[indx])
lsmI = lsfit(x=WI, y=WY, intercept=FALSE)
larm = lars(x=WX, y=WY, intercept=FALSE)
larm2 = lars(x=WX.2, y=WY, intercept=FALSE)
larm.c = lars(x=WX.c, y=WY.c, intercept=FALSE)

coef(larm)[,2:6] %*% cm
coef(larm)[,2:6] %*% wcm


XX = cbind(1,X)
as.matrix(coef(larm))->cl
XX %*% t(cl) -> ff 
apply(ff, 2, function(x) {sum(w[indx]*(Y-x)) / sum(w[indx])})