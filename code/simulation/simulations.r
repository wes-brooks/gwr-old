#library(sp, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(maps, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(shapefiles, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(plotrix, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(fossil, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(ggplot2, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))

#library(foreach, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(iterators, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(multicore, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(doMC, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#
#library(lars, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(glmnet, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(gwselect, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#
#library(splancs, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(RandomFields, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
#library(geoR, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))

library(geoR)
library(gwselect)
library(doMC)
registerDoMC(cores=7)

seeds = read.csv("seeds.csv")$x
N = 30
phi = 0.1
corr = 0.5
sigma = 1


coord = seq(0, 1, length.out=N)

args = commandArgs(trailingOnly=TRUE)
cluster = as.integer(args[1])
process = as.integer(args[2]) + 1

cluster=0
process=2

#Get two (independent) Gaussian random fields:
set.seed(seeds[process])
d1 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,phi))
d2 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,phi))
d3 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,phi))
d4 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,phi))
d5 = grf(n=N**2, grid='reg', cov.model='exponential', cov.pars=c(1,phi))

loc.x = d1$coords[,1]
loc.y = d1$coords[,2]

#Use the Cholesky decomposition to correlate the random fields:
S = matrix(corr, 5, 5)
diag(S) = rep(1, 5)
L = chol(S)

#Force correlation on the Gaussian random fields:
D = as.matrix(cbind(d1$data, d2$data, d3$data, d4$data, d5$data)) %*% L
    
#
X1 = matrix(D[,1], N, N)
X2 = matrix(D[,2], N, N)
X3 = matrix(D[,3], N, N)
X4 = matrix(D[,4], N, N)
X5 = matrix(D[,5], N, N)

B1 = matrix(rep( ifelse(coord<=0.4, 0, ifelse(coord<=0.6, 5*coord-2, 1)), N), N, N)
#B3 = matrix(rep(1-coord, N), N, N)



eta = X1*B1
#p = exp(eta) / (1+exp(eta))
Y = rnorm(N**2, eta, sigma)

#
loc.x = rep(seq(0, 1, length.out=N), each=N)
loc.y = rep(seq(0, 1, length.out=N), times=N)
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), X5=as.vector(X5), loc.x, loc.y)



#Find the optimal bandwidth and use it to generate a model:   
bw = gwlars.sel(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], mode.select="AIC", shrink=TRUE, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=TRUE, longlat=FALSE, adapt=TRUE, precondition=FALSE)
model = gwlars(Y~X1+X2+X3+X4+X5-1, data=sim, coords=sim[,c('loc.x','loc.y')], N=1, mode.select='AIC', shrink=TRUE, bw=bw, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=TRUE, longlat=FALSE, adapt=TRUE, precondition=FALSE)




#Write the results to some files:
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
for (k in 1:6) {
    coefs = t(sapply(1:N**2, function(y) {sapply(model[['model']][['models']][[y]][['coeflist']], function(x) {x[k]})}))
    write.table(coefs, file=paste("output/", vars[k], ".", cluster, ".", process, ".bootstrap.csv", sep=""), sep=',', row.names=FALSE, col.names=FALSE)
}

coefs = t(sapply(1:N**2, function(y) {as.vector(model[['model']][['models']][[y]][['coef']])}))
write.table(coefs, file=paste("output/CoefEstimates.", cluster, ".", process, ".csv", sep=""), col.names=vars, sep=',', row.names=FALSE)

params = c('bw', 'sigma2', 'loss.local', 's')
target = params[1]
output = sapply(1:N**2, function(y) {model[['model']][['models']][[y]][[target]]})

for (i in 2:length(params)) {
    target = params[i]
    output = cbind(output, sapply(1:N**2, function(y) {model[['model']][['models']][[y]][[target]]}))
}
write.table(output, file=paste("output/MiscParams.", cluster, ".", process, ".csv", sep=""), col.names=params, sep=',', row.names=FALSE)

