library(sp, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(shapefiles, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(plotrix, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(ggplot2, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(RandomFields, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(scales, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))

library(foreach, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(iterators, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(multicore, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(doMC, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))

library(lars, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(glmnet, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(gwselect, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))

library(splancs, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(geoR, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(maptools, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))
library(spgwr, lib.loc=c('R', 'R-libs/x86_64-redhat-linux-gnu-library/3.0'))


#library(geoR)
#library(gwselect)
#library(doMC)
#library(spgwr)
#registerCores(n=3)

seeds = as.vector(read.csv("seeds.csv", header=FALSE)[,1])
B = 100
N = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
settings = 6
tau = rep(0, settings)
rho = rep(c(0,0.5), settings/2)

b = 25

params = data.frame(tau, rho)

#Read command-line parameters
args = commandArgs(trailingOnly=TRUE)
cluster = as.integer(args[1])
process = as.integer(args[2])
#cluster=NA
#process=1

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


eta = X1*B1
Y = rbinom(n=length(eta), size=10, p=exp(eta)/(1+exp(eta)))/10

sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)
fitloc = cbind(rep(seq(0,1, length.out=N), each=N), rep(seq(0,1, length.out=N), times=N))

vars = cbind(B1=as.vector(B1!=0))#, B2=as.vector(B2!=0), B3=as.vector(B3!=0))
oracle = list()
for (i in 1:N**2) { 
    oracle[[i]] = character(0)
    if (vars[i,'B1']) { oracle[[i]] = c(oracle[[i]] , "X1") }
}


#Find the optimal bandwidth and use it to generate a model:
bw.glmnet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="BIC", gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, AICc=TRUE)
model.glmnet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='BIC', bw=bw.glmnet, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, AICc=TRUE)

bw.enet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="BIC", gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, AICc=TRUE)
model.enet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='BIC', bw=bw.enet, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, AICc=TRUE)

bw.oracular = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="BIC", gweight=bisquare, tol=0.01, method='dist', parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, AICc=TRUE)
model.oracular = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='binomial', oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='BIC', bw=bw.oracular, gweight=bisquare, tol=0.01, method='dist', simulation=TRUE, parallel=FALSE, interact=TRUE, verbose=FALSE, shrunk.fit=FALSE, AICc=TRUE)

bw.spgwr = ggwr.sel(Y~X1+X2+X3+X4+X5, family='binomial', data=sim, coords=as.matrix(sim[,c('loc.x','loc.y')]), gweight=gwr.bisquare, RMSE=FALSE, verbose=FALSE)
model.spgwr = ggwr(Y~X1+X2+X3+X4+X5, data=sim, family='binomial', coords=as.matrix(sim[,c('loc.x','loc.y')]), bandwidth=bw.spgwr, gweight=gwr.bisquare)




#First, write the data
write.table(sim, file=paste("output/Data.", cluster, ".", process, ".csv", sep=""), sep=',', row.names=FALSE)



#glmnet:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.glmnet[['model']][['models']][[y]][['coeflist']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".glmnet.bootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.glmnet[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".glmnet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


#Write the results to some files:
#vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.glmnet[['model']][['models']][[y]][['coef.unshrunk.list']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".glmnet.unshrunk-bootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.glmnet[['model']][['models']][[y]][['coef.unshrunk.interacted']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".unshrunk.glmnet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#ses = t(sapply(1:N**2, function(y) {as.vector(model.glmnet[['model']][['models']][[y]][['se.unshrunk']])}))
#write.table(ses, file=paste("output/CoefSEs.", cluster, ".", process, ".unshrunk.glmnet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


params = c('bw', 'sigma2', 'loss.local', 's', 's2.unshrunk', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.glmnet[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.glmnet[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".glmnet.csv", sep=""), col.names=params, sep=',', row.names=FALSE)



#enet:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.enet[['model']][['models']][[y]][['coeflist']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".enet.bootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


#Write the results to some files:
#vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.enet[['model']][['models']][[y]][['coef.unshrunk.list']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".enet.unshrunk-bootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['coef.unshrunk.interacted']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".unshrunk.enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#ses = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['se.unshrunk']])}))
#write.table(ses, file=paste("output/CoefSEs.", cluster, ".", process, ".unshrunk.enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


params = c('bw', 'sigma2', 'loss.local', 's', 's2.unshrunk', 'fitted')
target = params[1]
output = sapply(model.enet[['model']][['models']], function(y) {y[[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(model.enet[['model']][['models']], function(y) {y[[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".enet.csv", sep=""), col.names=params, sep=',', row.names=FALSE)








#For oracle property:
#Write the results to some files:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.oracular[['model']][['models']][[y]][['coeflist']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".OracularBootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.oracular[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".oracle.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#ses = t(sapply(1:N**2, function(y) {as.vector(model.oracular[['model']][['models']][[y]][['se.coef']])}))
#write.table(ses, file=paste("output/CoefSEs.", cluster, ".", process, ".oracle.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


params = c('bw', 'sigma2', 'loss.local', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.oracular[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.oracular[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".oracle.csv", sep=""), col.names=params, sep=',', row.names=FALSE)









#For all vars (spgwr):
#Write the results to some files:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = as.matrix(model.spgwr$SDF@data[,2:7])
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".spgwr.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

d = cbind(1, sim[,2:6])
l = diag(as.matrix(d) %*% t(as.matrix(model.spgwr$SDF@data[,2:7])))
f = exp(l) / (1+exp(l))
output = as.matrix(f)
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".spgwr.csv", sep=""), col.names=c("fitted"), sep=',', row.names=FALSE)
