library(gwselect)
library(MASS)
library(geoR)
registerCores(n=7)

size = c(20, 30, 40, 50)
size = c(40)

for (N in size) {
    coord = seq(0, 1, length.out=N)
    grid = matrix(rnorm(N**2, mean=rep(ifelse(coord<=0.5, 0, 2), N), sd=0.5), N, N)
    rownames(grid) = seq(0, 1, length.out=N)
    colnames(grid) = seq(0, 1, length.out=N)
}

#population for weights
pop = rpois(N**2, 400)

#
d1 = mvrnorm(n=N**2, mu=c(0,0), Sigma=matrix(c(1,0.2,0.2,1),2,2))
d2 = mvrnorm(n=N**2, mu=c(0,0), Sigma=matrix(c(1,0.2,0.2,1),2,2))

x1 = d1[,1]
X1 = matrix(x1, N, N)
B1 = matrix(rep(ifelse(coord<=0.5, 0, 2), N), N, N)

x2 = d2[,1]
X2 = matrix(x2, N, N)
B2 = matrix(rep(1-coord, N), N, N)

#Correlated with X1:
x3 = d1[,2]
X3 = matrix(x3, N, N)

#Correlated with X2
x4 = d2[,2]
X4 = matrix(x4, N, N)

eta = X1*B1 + X2*B2 + rnorm(N**2, 0, 0.1)
Z = rnorm(N**2, 0, 1)
#p = exp(eta) / (1+exp(eta))
Y = rnorm(N**2, eta, 1)

#
loc.x = rep(seq(0, 1, length.out=N), each=N)
loc.y = rep(seq(0, 1, length.out=N), times=N)
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), X3=as.vector(X3), X4=as.vector(X4), Z, loc.x, loc.y)

#weights = pop
#bw = gwglmnet.sel(Y~X1+X2+X3+X4+Z, data=sim, coords=sim[,c('loc.x','loc.y')], weights=pop, gweight=bisquare, tol=0.01, s=NULL, method='knn', family='binomial', parallel=TRUE, longlat=FALSE, adapt=TRUE, precondition=FALSE)
#model = gwglmnet(Y~X1+X2+X3+X4+Z, data=sim, coords=sim[,c('loc.x','loc.y')], bw=bw, weights=pop, gweight=bisquare, tol=0.01, s=NULL, method='knn', family='binomial', parallel=FALSE, longlat=FALSE, adapt=TRUE, precondition=FALSE)

bw = gwlars.sel(Y~X1+X2+X3+X4+Z, data=sim, coords=sim[,c('loc.x','loc.y')], mode.select="AIC", weights=rep(1, N**2), gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=TRUE, longlat=FALSE, adapt=TRUE, precondition=FALSE)
#Bandwidth: 0.00636730579145969. CV error: 966.753380263565
bw.precondition = gwlars.sel(Y~X1+X2+X3+X4+Z, data=sim, coords=sim[,c('loc.x','loc.y')], weights=pop, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=TRUE, longlat=FALSE, adapt=FALSE, precondition=TRUE)
model = gwlars(Y~X1+X2+X3+X4+Z, data=sim, coords=sim[,c('loc.x','loc.y')], mode.select='AIC', bw=bw, weights=pop, gweight=bisquare, tol=0.01, s=NULL, mode='step', method='knn', parallel=FALSE, longlat=FALSE, adapt=TRUE, precondition=FALSE)



#model = gwglmnet.nen(nifestations~meanelevation+warm+Tmin+Tmean+Tmax+cold+precip+dd+ddegg, data=mpb, coords=mpb[,c('X','Y')], gweight=bisquare, s=seq(0,5,0.001), tol=10, bw=200000, type='pearson', family='binomial', parallel=TRUE, weights=weights)
#0.004797285
#0.00479728530129436