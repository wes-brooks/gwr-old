#Integrating the product of two thin-plate spline basis functions:
library(spgwr)
library(lars)
library(mgcv)
library(spgwr)
data(columbus)


bisquare = function(zi, zj, bandwidth) {
    ifelse(r(zi[1], zi[2], zj) < bandwidth, (1 - (r(zi[1], zi[2], zj)/bandwidth)**2)**2, 0)
}


r = function(x, y, zk) { 
    sqrt((x-zk[1])**2 + (y-zk[2])**2)
}

#Hardcoded for two dimensions:
eta = function(x, y, zk) {
    ifelse(x==zk[1] & y==zk[2], 0, r(x, y, zk)**2 * log(r(x, y, zk)) / (8 * pi) )
}


f = function(y, x, zi, zj) {
    eta(x, y, zi) * eta(x, y, zj)
}


fvec = function(y, x, zi, zj) {
    sapply(y, f, x=x, zi=zi, zj=zj) 
}  


g = function(x, zi, zj, ylim) {
    integrate(fvec, lower=ylim[1], upper=ylim[2], x=x, zi=zi, zj=zj)$val
}


gvec = function(x, zi, zj, ylim) {
    sapply(x, g, zi=zi, zj=zj, ylim=ylim)
} 

zi = c(0.2, 0.2)
zj = c(0.5, 0.5)

xlim = range(columbus$x)
ylim = range(columbus$y)

#integrate(gvec, lower=xlim[1], upper=xlim[2], zi=zi, zj=zj, ylim=ylim) 

#m1 = lm(housing ~ income + crime, data=columbus)
#columbus$resid = m1$residuals

crime = vector()
income = vector()
intercept = vector()

int.l = vector()
inc.l = vector()
cr.l = vector()
w.lasso = list()

for(i in 1:dim(columbus)[1]) {
    w = vector()
    for(j in 1:dim(columbus)[1]) {
        w = c(w, bisquare(c(columbus$x[i], columbus$y[i]), c(columbus$x[j], columbus$y[j]), bandwidth=5)) }

    model = lm(housing ~ income + crime, data=columbus, weights=w)
    
    w.eig <- eigen(diag(w))
    w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    w.lasso[[i]] = lars(x=w.sqrt %*% as.matrix(columbus[,c("crime", "income")]), y=as.matrix(columbus[,"housing"]))

    intercept = c(intercept, model$coef[["(Intercept)"]])
    income = c(income, model$coef[["income"]])
    crime = c(crime, model$coef[["crime"]])
}

x = columbus$x
y = columbus$y

#Create a thin plate spline that interpolates the GWR model's coefficients:
tps = gam(income ~ s(x, y, k=49))


E = matrix(0, 49, 49)
for (i in 1:49) {
    for (j in 1:49) {
        E[i,j] = eta(x[i], y[i], c(x[j],y[j]))
    }
}

T = cbind(1, x, y)
U = eigen(E)[["vectors"]][,1:48]
TU = t(T) %*% U



#Find the matrix of the integrated products of thin plate spline basis functions:
R = matrix(NA, nrow=dim(columbus)[1], ncol=dim(columbus)[1])
for(i in 1:dim(columbus)[1]) {

    zi = c(columbus$x[i], columbus$y[i])

    for(j in i:dim(columbus)[1]) {

        zj = c(columbus$x[j], columbus$y[j])
        R[i,j] = integrate(gvec, lower=xlim[1], upper=xlim[2], zi=zi, zj=zj, ylim=ylim)$val
    }
}

for(i in 1:(dim(columbus)[1]-1)) {
    for(j in (i+1):dim(columbus)[1]) {
        R[j,i] = R[i,j]
    }
}
