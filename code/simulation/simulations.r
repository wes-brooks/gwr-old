library(sp, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(maps, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(shapefiles, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(plotrix, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(fossil, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(ggplot2, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))

library(foreach, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(iterators, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(multicore, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(doMC, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))

library(lars, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(glmnet, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(gwselect, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))

library(splancs, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(RandomFields, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
library(geoR, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))

#library(geoR)
#library(gwselect)
#library(doMC)
#registerDoMC(cores=7)

seeds = read.csv("seeds.csv")$x
#size = c(20, 30, 40, 50)
N = 30
coord = seq(0, 1, length.out=N)

args = commandArgs(trailingOnly=TRUE)
cluster = args[1]
process = args[2]

#Get two (independent) Gaussian random fields:
set.seed(seeds[process])
d1a = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,1))
d1b = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,1))
loc.x = d1a$coords[,1]
loc.y = d1a$coords[,2]

#Use the Cholesky decomposition to correlate the two fields:
S1 = matrix(c(1,0.2,0.2,1),2,2)
L1 = chol(S1)

#Force correlation on the Gaussian random fields:
D1 = as.matrix(cbind(d1a$data, d1b$data)) %*% L1
    
#Get two (independent) Gaussian random fields:
d2a = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,1))
d2b = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,1))  

#Use the Cholesky decomposition to correlate the two fields:
S2 = matrix(c(1,0.2,0.2,1),2,2)
L2 = chol(S2)

#Force correlation on the Gaussian random fields:
D2 = as.matrix(cbind(d2a$data, d2b$data)) %*% L2
    
#
X1 = matrix(D1[,1], N, N)
B1 = matrix(rep(ifelse(coord<=0.5, 0, 2), N), N, N)

X2 = matrix(D1[,2], N, N)
    
#
X3 = matrix(D2[,1], N, N)
B3 = matrix(rep(1-coord, N), N, N)

#Correlated with X3
X4 = matrix(D2[,2], N, N)
    
d3 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,1)) 
X5 = matrix(d3$data, N, N)
#B5 = t(matrix(rep(3-3*coord, N), N, N))

eta = X1*B1 + X3*B3 #+ X5*B5
Z = rnorm(N**2, 0, 1)
#p = exp(eta) / (1+exp(eta))
Y = rnorm(N**2, eta, 1)

#
loc.x = rep(seq(0, 1, length.out=N), each=N)
loc.y = rep(seq(0, 1, length.out=N), times=N)
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), Z, loc.x, loc.y)

#Find the optimal bandwidth and use it to generate a model:   
bw = gwlars.sel(Y~X1+X2+X3+X4+X5+Z, data=sim, coords=sim[,c('loc.x','loc.y')], mode.select="AIC", gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=FALSE, longlat=FALSE, adapt=TRUE, precondition=FALSE)
model = gwlars(Y~X1+X2+X3+X4+X5+Z-1, data=sim, coords=sim[,c('loc.x','loc.y')], N=101, mode.select='AIC', bw=bw, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=FALSE, longlat=FALSE, adapt=TRUE, precondition=FALSE)

#Write the results to some files:
vars = c('Intercept', 'X1', 'X2', 'X3', 'X4', 'X5', 'Z')
for (k in 2:7) {
    coefs = t(sapply(1:N**2, function(y) {sapply(model[['model']][['models']][[y]][['coeflist']], function(x) {x[k]})}))
    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".bootstrap.csv", sep=""))
}

coefs = t(sapply(1:N**2, function(y) {sapply(model[['model']][['models']][[y]][['intlist']], function(x) {x[1]})}))
write.table(coefs, file=paste("output/", Intercept, ".", cluster, ".", process, ".bootstrap.csv", sep=""))

coefs = t(sapply(1:N**2, function(y) {as.vector(model[['model']][['models']][[y]][['coef']])}))
coefs = cbind(coefs, sapply(1:N**2, function(y) {model[['model']][['models']][[y]][['intercept']]}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".csv", sep=""), col.names=vars)


params = c('bw', 'sigma2', 'loss.local', 's', 'sum.weights')
target = params[1]
output = sapply(1:N**2, function(y) {model[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".csv", sep=""), col.names=params)



#Bandwidth: 0.124453494228221
#Bandwidth: 0.00636730579145969. CV error: 966.753380263565
#bw=0.0502464695551903
#bw=0.05572809
#bw.precondition = gwlars.sel(Y~X1+X2+X3+X4+Z, data=sim, coords=sim[,c('loc.x','loc.y')], weights=pop, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=TRUE, longlat=FALSE, adapt=FALSE, precondition=TRUE)
#model = gwlars(Y~X1+X2+X3+X4+X5+Z-1, data=sim, coords=sim[,c('loc.x','loc.y')], N=101, mode.select='AIC', bw=bw, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=TRUE, longlat=FALSE, adapt=TRUE, precondition=FALSE)

#Weights set to 1, weighting moved to front of algorithm
#Bandwidth: 0.00402385017589017

#model = gwglmnet.nen(nifestations~meanelevation+warm+Tmin+Tmean+Tmax+cold+precip+dd+ddegg, data=mpb, coords=mpb[,c('X','Y')], gweight=bisquare, s=seq(0,5,0.001), tol=10, bw=200000, type='pearson', family='binomial', parallel=TRUE, weights=weights)
#0.004797285
#0.00479728530129436