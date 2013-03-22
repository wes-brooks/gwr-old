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
