require(survival)
data(kidney)
kidney$frail = 10*kidney$frail

bw = 



    #Isolate one year of data
    year = as.character(yr)
    df = pov2[pov2$year==yr,]
    
    ####Produce the models via lasso and via elastic net:
    #Define which variables we'll use as predictors of poverty:
    #predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')
    f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))
    bw[['GWAL']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, mode.select="BIC", gweight=spherical, tol=0.1, s=NULL, method='dist', adapt=TRUE, parallel=FALSE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)
    model[['GWAL']][[year]] = gwglmnet(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, N=1, mode.select='BIC', bw=bw[['GWAL']][[year]][['res']], gweight=spherical, method='dist', simulation=TRUE, adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)


trace = bw[['GWAL']][[year]][['trace']]
bw.mle = bw[['GWAL']][[year]][['res']]

a = diag(rep(-1,22))
aa = diag(rep(1,21))
a[2:22,1:21] = a[2:22,1:21] + aa

(t(trace[,1]) %*% a)[1:21] -> run 
(t(trace[,2]) %*% a)[1:21] -> rise 


b = diag(rep(-1,21))
bb = diag(rep(1,20))
b[2:21,1:20] = b[2:21,1:20] + bb

(rise %*% b)[1:20] -> rise2 
(run %*% b)[1:20] -> run2 

#A weighted mean:
c1 = mean(-trace[1:20,1]*rise2/run2/sum(trace[1:20,1]))

#Or just grab from the region around the MLE
c2 = 0.0272

xx = c(1300,1500)
yy = c(-665,-655)
plot(trace[4:22,], xlim=xx, ylim=yy, bty='n')

xxx = seq(xx[1],xx[2], length.out=200)
yyy1 = c1*(xxx - bw.mle)**2 + min(trace[,2])
par(new=TRUE)
plot(xxx,yyy1, type='l', lty=2, col='red', xlim=xx, ylim=yy, xaxt='n', yaxt='n', ann=FALSE, bty='n')

yyy2 = c2*(xxx - bw.mle)**2 + min(trace[,2])
par(new=TRUE)
plot(xxx,yyy2, type='l', lty=3, col='green', xlim=xx, ylim=yy, xaxt='n', yaxt='n', ann=FALSE, bty='n')



converting to a gaussian density:
