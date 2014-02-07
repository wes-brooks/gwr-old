library(sp, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(shapefiles, lib.loc'R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(plotrix, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(ggplot2, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(RandomFields, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(scales, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')

library(foreach, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(iterators, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(multicore, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(doMC, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')

library(lars, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(glmnet, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(gwselect, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')

library(splancs, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(geoR, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(maptools, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')
library(spgwr, lib.loc='R-libs/x86_64-redhat-linux-gnu-library/3.0')


seeds = as.vector(read.csv("seeds.csv", header=FALSE)[,1])
B = 100
N = 30
settings = 9
functions = 3
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(0.1, settings)
rho = rep(c(0, 0.5, 0.9), functions)
params = data.frame(tau, rho)

#Read command-line parameters
args = commandArgs(trailingOnly=TRUE)
cluster = as.integer(args[1])
process = as.integer(args[2])

#Simulation parameters are based on the value of process
setting = process %/% B + 1
parameters = params[setting,]
set.seed(seeds[process+1])

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

if ((setting-1) %/% 3 == 0) {
    B1 = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
} else if ((setting-1) %/% 3 == 1) {
    B1 = matrix(rep(coord, N), N, N)
} else if ((setting-1) %/% 3 == 2) {
    Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
    Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
    D = (Xmat-0.5)**2 + (Ymat-0.5)**2
    d = D[,435]
    B1 = matrix(max(d)-d, N, N)
}

m = 1
eta = X1*B1
Y = rbinom(n=length(eta), size=m, p=exp(eta)/(1+exp(eta))) / m

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
fitloc = cbind(rep(seq(0,1, length.out=N), each=N), rep(seq(0,1, length.out=N), times=N))

vars = cbind(B1=as.vector(B1!=0))#, B2=as.vector(B2!=0), B3=as.vector(B3!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}


#MODELS:

bw.glmnet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select='AICc', gweight=spherical, tol.bw=0.01, bw.method='knn', precondition=FALSE, parallel=FALSE, interact=FALSE, verbose=TRUE, shrunk.fit=TRUE, bw.select='AICc', resid.type='pearson')
model.glmnet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AICc', bw=bw.glmnet[['bw']], gweight=spherical, bw.method='knn', simulation=TRUE, precondition=FALSE, parallel=FALSE, interact=FALSE, verbose=TRUE, shrunk.fit=TRUE)

bw.enet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select='AICc', gweight=spherical, tol.bw=0.01, bw.method='knn', precondition=FALSE, parallel=FALSE, interact=FALSE, verbose=TRUE, shrunk.fit=TRUE, bw.select='AICc', resid.type='pearson')
model.enet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AICc', bw=bw.enet[['bw']], gweight=spherical, bw.method='knn', simulation=TRUE, precondition=FALSE, parallel=FALSE, interact=FALSE, verbose=TRUE, shrunk.fit=TRUE)

bw.oracular = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, gweight=spherical, tol.bw=0.01, bw.method='knn', parallel=FALSE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE, bw.select='AICc', resid.type='pearson')
model.oracular = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, bw=bw.oracular[['bw']], gweight=spherical, bw.method='knn', simulation=TRUE, parallel=FALSE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

oracle2 = lapply(1:900, function(x) {return(c("X1", "X2", "X3", "X4", "X5"))})
bw.gwr = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', oracle=oracle2, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, gweight=spherical, tol.bw=0.01, bw.method='knn', parallel=FALSE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE, bw.select='AICc', resid.type='pearson')
model.gwr = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', oracle=oracle2, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, bw=bw.oracular[['bw']], gweight=spherical, bw.method='knn', simulation=TRUE, parallel=FALSE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)


#OUTPUT:

#First, write the data
write.table(sim, file=paste("output/Data.", cluster, ".", process, ".csv", sep=""), sep=',', row.names=FALSE)


#glmnet:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = t(sapply(1:N**2, function(y) {as.vector(model.glmnet[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".glmnet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

coefs = t(sapply(1:N**2, function(y) {as.vector(model.glmnet[['model']][['models']][[y]][['coef.unshrunk']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".unshrunk.glmnet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

params = c('bw', 'sigma2', 'loss.local', 's2.unshrunk', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.glmnet[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.glmnet[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".glmnet.csv", sep=""), col.names=params, sep=',', row.names=FALSE)






#enet:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

coefs = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['coef.unshrunk']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".unshrunk.enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

params = c('bw', 'sigma2', 'loss.local', 's2.unshrunk', 'fitted')
target = params[1]
output = sapply(model.enet[['model']][['models']], function(y) {y[[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(model.enet[['model']][['models']], function(y) {y[[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".enet.csv", sep=""), col.names=params, sep=',', row.names=FALSE)








#For oracle property:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = t(sapply(1:N**2, function(y) {as.vector(model.oracular[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".oracle.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

params = c('bw', 'sigma2', 'loss.local', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.oracular[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.oracular[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".oracle.csv", sep=""), col.names=params, sep=',', row.names=FALSE)





#For all vars:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = t(sapply(1:N**2, function(y) {as.vector(model.gwr[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".gwr.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

params = c('bw', 'sigma2', 'loss.local', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.gwr[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.gwr[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".gwr.csv", sep=""), col.names=params, sep=',', row.names=FALSE)

