library(sp, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(shapefiles, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(plotrix, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(ggplot2, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(RandomFields, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(scales, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))

library(foreach, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(iterators, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(multicore, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(doMC, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))

library(lars, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(glmnet, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(gwselect, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))

library(splancs, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(geoR, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(maptools, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))
library(spgwr, lib.loc=c('R', 'R-libs/i386-redhat-linux-gnu-library'))

#library(geoR)
#library(gwselect)
#library(doMC)
#registerCores(n=3)

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
args = commandArgs(trailingOnly=TRUE)
cluster = as.integer(args[1])
process = as.integer(args[2])
#cluster=NA
#process=2

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
bw.lars = gwlars.sel(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,1), gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)
model.lars = gwlars(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.lars, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)

bw.glmnet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,1), gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)
model.glmnet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha=1, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.glmnet, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)

bw.enet = gwglmnet.sel(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,1), gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)
model.enet = gwglmnet(Y~X1+X2+X3+X4+X5-1, data=sim, family='gaussian', alpha='adaptive', coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.enet, gweight=bisquare, tol=0.01, s=NULL, method='dist', simulation=TRUE, adapt=TRUE, precondition=FALSE, parallel=FALSE, interact=TRUE)

bw.oracular = gwlars.sel(Y~X1+X2+X3+X4+X5-1, data=sim, oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,1), gweight=bisquare, tol=0.01, method='dist', parallel=FALSE, interact=TRUE)
model.oracular = gwlars(Y~X1+X2+X3+X4+X5-1, data=sim, oracle=oracle, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.oracular, gweight=bisquare, tol=0.01, method='dist', simulation=TRUE, parallel=FALSE, interact=TRUE)

oracle2 = lapply(1:900, function(x) {return(c("X1", "X2", "X3", "X4", "X5"))})
bw.gwr = gwlars.sel(Y~X1+X2+X3+X4+X5-1, data=sim, oracle=oracle2, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, mode.select="AIC", range=c(0,1), gweight=bisquare, tol=0.01, method='dist', parallel=FALSE, interact=FALSE)
model.gwr = gwlars(Y~X1+X2+X3+X4+X5-1, data=sim, oracle=oracle2, coords=sim[,c('loc.x','loc.y')], longlat=FALSE, N=1, mode.select='AIC', bw=bw.gwr, gweight=bisquare, tol=0.01, method='dist', simulation=TRUE, parallel=FALSE, interact=FALSE)

bw.spgwr = gwr.sel(Y~X1+X2+X3+X4+X5, data=sim, coords=as.matrix(sim[,c('loc.x','loc.y')]), gweight=gwr.bisquare, method="aic")
model.spgwr = gwr(Y~X1+X2+X3+X4+X5, data=sim, coords=as.matrix(sim[,c('loc.x','loc.y')]), bandwidth=bw.spgwr, gweight=gwr.bisquare)

#First, write the data
write.table(sim, file=paste("output/Data.", cluster, ".", process, ".csv", sep=""), sep=',', row.names=FALSE)




#LARS:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.lars[['model']][['models']][[y]][['coeflist']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".bootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.lars[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".lars.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


#Write the results to some files:
#vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
#for (k in 1:6) {
#    coefs = t(sapply(1:N**2, function(y) {sapply(model.lars[['model']][['models']][[y]][['coef.unshrunk.list']], function(x) {x[k]})}))
#    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".unshrunk-bootstrap.csv", sep=""), sep=',', row.names=FALSE)
#}

coefs = t(sapply(1:N**2, function(y) {as.vector(model.lars[['model']][['models']][[y]][['coef.unshrunk']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".unshrunk.lars.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#ses = t(sapply(1:N**2, function(y) {as.vector(model.lars[['model']][['models']][[y]][['se.unshrunk']])}))
#write.table(ses, file=paste("output/CoefSEs.", cluster, ".", process, ".unshrunk.lars.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


params = c('bw', 'sigma2', 'loss.local', 's', 's2.unshrunk', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.lars[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.lars[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".lars.csv", sep=""), col.names=params, sep=',', row.names=FALSE)





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

coefs = t(sapply(1:N**2, function(y) {as.vector(model.glmnet[['model']][['models']][[y]][['coef.unshrunk']])}))
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

coefs = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['coef.unshrunk']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".unshrunk.enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

#ses = t(sapply(1:N**2, function(y) {as.vector(model.enet[['model']][['models']][[y]][['se.unshrunk']])}))
#write.table(ses, file=paste("output/CoefSEs.", cluster, ".", process, ".unshrunk.enet.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


params = c('bw', 'sigma2', 'loss.local', 's', 's2.unshrunk', 'fitted')
target = params[1]
output = sapply(1:N**2, function(y) {model.enet[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model.enet[['model']][['models']][[y]][[target]]}))
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






#For all vars (gwlars.oracle):
#Write the results to some files:
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







#For all vars (spgwr):
#Write the results to some files:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

coefs = as.matrix(model.spgwr$SDF@data[,2:7])
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".spgwr.csv", sep=""), col.names=vars, sep=',', row.names=FALSE)


output = as.matrix(model.spgwr$SDF@data[,'pred'])
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".spgwr.csv", sep=""), col.names=c("fitted"), sep=',', row.names=FALSE)

