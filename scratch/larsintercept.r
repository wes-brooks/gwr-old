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
XX = cbind(1,X)
Y = as.matrix(sim$Y[indx])
wm = sum(w[indx]*Y) / sum(w)

WX = W %*% X
WXX = W%*%XX
WY = W %*% Y

WY.c = WY-mean(WY)
WX.c = scale(WX, scale=FALSE)
WX.cs = scale(WX)
mean.wx = colMeans(WX)
norm.wx = attr(WX.cs, "scaled:scale")

WX.c = cbind(sqrt(w[indx]), WX.c)
WX.cs = cbind(sqrt(w[indx]), WX.cs)

cm = apply(X, 2, function(x) {sum(x*w[indx])})/sum(w)
wcm = apply(WX, 2, function(x) {sum(x*sqrt(w[indx]))})/sum(sqrt(w))

lsm = lsfit(x=X, y=Y, wt=w[indx])
larm = lars(x=WXX, y=WY, intercept=FALSE, normalize=FALSE)
larm.c = lars(x=WX.c, y=WY, intercept=FALSE, normalize=FALSE)
larm.cc = lars(x=WX.c, y=WY.c, intercept=FALSE, normalize=FALSE)
larm.cs = lars(x=WX.cs, y=WY, intercept=FALSE, normalize=FALSE)

coef(larm)[,2:6] %*% cm
coef(larm)[,2:6] %*% wcm



as.matrix(coef(larm))->cl
XX %*% t(cl) -> ff 
apply(ff, 2, function(x) {sum(w[indx]*(Y-x)) / sum(w[indx])})


as.matrix(coef(larm.cs))->cl.cs
cl.cs.rescaled = t(apply(cl.cs, 1, function(x) {x * c(1, 1/norm.wx)}))

coef(larm.c)[7,1]-sum(coef(larm.c)[7,2:6] * cm)
apply(coef(larm.c), 1, function(x) {x[1] - sum(x[-1] * cm)})

as.matrix(coef(larm.cs))->cl.cs
XX.c = 
XX %*% t(cl.cs) -> ff 
apply(ff, 2, function(x) {sum(w[indx]*(Y-x)) / sum(w[indx])})
