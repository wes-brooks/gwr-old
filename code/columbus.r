library(spgwr)
data(columbus)

#Import the plotting functions:
setwd("~/git/gwr/code")
source("utils.r")

#Define a grid of locations where we'll fit a GWR model:
n=20
xx = as.vector(quantile(columbus$x, 1:n/(n+1)))
yy = as.vector(quantile(columbus$y, 1:n/(n+1)))
locs = cbind(x=rep(xx,each=n), y=rep(yy,times=n))

#Use this trick to compute the matrix of distances very quickly
n = dim(columbus)[1]
D1 = matrix(rep(columbus$x,n), n,n)
D2 = matrix(rep(columbus$y,n), n,n)
D = sqrt((D1-t(D1))**2 + (D2-t(D2))**2)

#Define which variables we'll use as predictors of house prices:
predictors = c('crime', 'income')
f = as.formula(paste("housing ~ ", paste(predictors, collapse="+"), sep=""))

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}

#Use the lasso for GWR models of income with 2006 data:
df = columbus
w.lasso.geo = list()
coefs = list()

for(i in 1:dim(df)[1]) {
    w = bisquare(D[,i], bw=bw)

    model = lm(f, data=df, weights=w)
    
    #w.eig <- eigen(diag(w))
    #w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    #w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[,predictors]), y=as.matrix(df$logitindpov))
    
    for (col in predictors) {
        coefs[[col]] = c(coefs[[col]], model$coef[[col]])
    }
    
    print(i)
}


#Use the methods of spgwr to select a bandwidth and fit a GWR model for poverty:
bw = gwr.sel(housing~crime+income, data=columbus, coords=cbind(x,y), adapt=FALSE, gweight=gwr.bisquare)
gwr.columbus = gwr(housing~crime+income, data=columbus, coords=cbind(x,y), bandwidth=bw, gweight=gwr.bisquare, hatmatrix=TRUE)
