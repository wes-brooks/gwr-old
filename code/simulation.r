#Simulate data for testing a GWR model:
setwd("~/git/gwr/code")
source("utils.r")
library(lars)
library(spgwr)

#Set up the location grid:
N = 20
x = rep(0:N - N/2, each=N+1)
y = rep(0:N - N/2, times=N+1)

#Simulate the predictors independently of location:
mu.A = 0
sig.A = 1
A = rnorm((N+1)**2, mu.A, sig.A)

mu.B = 1
sig.B = 2
B = rnorm((N+1)**2, mu.B, sig.B)

mu.C = -1
sig.C = 3
C = rnorm((N+1)**2, mu.C, sig.C)

#Simulate the output variable:
sig.err = 0.5
out = y + A*x + B + rnorm((N+1)**2, 0, sig.err)

simulated = data.frame(x, y, A, B, C, out)

#Use the methods of spgwr to select a bandwidth and fit a GWR model for poverty:
bw = gwr.sel(out~A+B+C, data=simulated, coords=cbind(x,y), adapt=FALSE, gweight=gwr.bisquare)
gwr.sim = gwr(out~A+B+C, data=simulated, coords=cbind(x,y), bandwidth=bw, gweight=gwr.bisquare)



#Homebrew GWR:
df = simulated

#Define which variables we'll use as predictors of poverty:
predictors = c('A', 'B', 'C')
output = 'out'
f = as.formula(paste("out ~ 1 + ", paste(predictors, collapse="+"), sep=""))

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}

#Use this trick to compute the matrix of distances very quickly
n = dim(df)[1]
D1 = matrix(rep(df$x,n), n,n)
D2 = matrix(rep(df$y,n), n,n)
D = sqrt((D1-t(D1))**2 + (D2-t(D2))**2)

#Use the lasso for GWR model selection:
w.lasso.geo = list()
coefs = list()
diagnostics = list()

for(i in 1:dim(df)[1]) {
    w = bisquare(D[,i], bw=bw)

    model = lm(f, data=df, weights=w)
    
    w.eig <- eigen(diag(w))
    w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[,predictors]), y=as.matrix(df[[output]]))
    
    for (col in predictors) {
        coefs[[col]] = c(coefs[[col]], model$coef[[col]])
    }
    coefs[['(Intercept)']] = c(coefs[['(Intercept)']], model$coef[['(Intercept)']])

    diagnostics[['R2']] = c(diagnostics[['R2']], summary(model)[['r.squared']])
    diagnostics[['sigma']] = c(diagnostics[['sigma']], summary(model)[['sigma']])
    diagnostics[['total weight']] = c(diagnostics[['total weight']], sum(w))
    
    print(i)
}

model = list(coords=df[,c('x','y')], coefs=coefs, diags=diagnostics, lasso=w.lasso.geo)

