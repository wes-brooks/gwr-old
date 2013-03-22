n=900
coords = sim[,c('loc.x', 'loc.y')]
Xmat = matrix(rep(coords[,1], times=n), n, n)
Ymat = matrix(rep(coords[,2], times=n), n, n)
D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
            
d = D[320,]
bw = 0.2
w = bisquare(d, bw)
indx = which(w>0)

X = as.matrix(sim[indx,c('X1','X2','X3','X4','X5')])
X2 = cbind(1,X)
Y = as.matrix(sim$Y[indx])

W = diag(sqrt(w[indx]))

WX = W%*%X
WY = W%*%Y
wm = sum(WY)/sum(W)
m = lars(x=WX, y=WY, type='lar')

WX2 = WX - colMeans(WX)
mm = lars(x=WX2, y=WY2, type='lar', normalize=FALSE)
mmm = lars(x=WX2, y=WY2, type='lar', normalize=FALSE, intercept=FALSE)

X3 = cbind(1, X)
WX3 = W%*%X3
WX3.c = W%*%X3
WX3.c[,2:ncol(WX3)] = WX3[,2:ncol(WX3)] - colMeans(WX3[,2:ncol(WX3)])

WY.c = WY - mean(WY)
m3 = lars(x=WX3.c, y=WY, normalize=FALSE, intercept=FALSE)


m1 = lsfit(x=X, y=Y, wt=w[indx])

mg = glmnet(x=X, y=Y, family='gaussian', weights=w[indx])