mpb = read.csv("~/git/gwr/data/SouthernPineBeetle/Code-Andy/mpb.csv", head=TRUE)

source("~/git/gwr/code/utils.r")
source("~/git/gwr/code/PineBeetle/nearest_effective_neighbors.r")


m = gwglmnet.nen(nifestations~meanelevation+warm, data=mpb, coords=mpb[,c('X','Y')], gweight=bisquare, s=seq(0,5,0.001), tol=10, type='pearson', bw=200000, family='poisson', parallel=TRUE)


#gwglmnet.sel(formula=nifestations~meanelevation+Tmin+Tmean+Tmax+warm+cold+precip+dd+ddegg, data=mpb, family='poisson', coords=mpb[,c('X','Y')], longlat=FALSE, adapt=TRUE, gweight=bisquare, s=seq(0,5,0.001), tol=0.1)

