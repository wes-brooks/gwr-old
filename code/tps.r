#Integrating the product of two thin-plate spline basis functions:
library(spgwr)
library(lars)
library(mgcv)
library(spgwr)
data(columbus)


bisquare = function(x, z, bandwidth) {
    ifelse( r(x=x, z=c(z[1],z[2])) < bandwidth, (1 - (r(x=x, z=c(z[1],z[2]))/bandwidth)**2)**2, 0)
}


r = function(x, z) { 
    sqrt((x[1]-z[1])**2 + (x[2]-z[2])**2)
}

#Hardcoded for two dimensions:
R = function(x, z) {
    ifelse(x[1]==z[1] & x[2]==z[2], 0, r(x, z)**2 * log(r(x, z)) / (8 * pi) )
}


f = function(y, x, zi, zj) {
    R(c(x, y), zi) * R(c(x, y), zj)
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






spl.X <- function(x) {
    q = length(unique(x)) + 2
    n = dim(x)[1]

    for(col in 1:dim(x)[2]) {
        x[,col] = x[,col] - mean(x[,col])
    }    

    E = matrix(0, n, n)
    for (i in 1:n) {
        for (j in 1:n) {
            xi = as.vector(c(x[i,1], x[i,2]))
            xj = as.vector(c(x[j,1], x[j,2]))
            E[i,j] = R(x=xi, z=xj)
        }
    }
    
    X = cbind(rep(1,n), x, E)
    return(X)
}




major = order(abs(eigen.decomp$values), decreasing=TRUE)[1:k]
deflated = eigen.decomp$vectors[,major] %*% diag(eigen.decomp$values[major]) %*% t(eigen.decomp$vectors[,major])



T = cbind(1, columbus$x-mean(columbus$x), columbus$y-mean(columbus$y))

n = dim(columbus)[1]
E = matrix(0, n, n)

for (i in 1:n) {
    for (j in 1:n) {
        xi = as.vector(c(columbus[i,'x'], columbus[i,'y']))
        xj = as.vector(c(columbus[j,'x'], columbus[j,'y']))
        E[i,j] = R(x=xi, z=xj)
    }
}

k=48
eigen.decomp = eigen(E)
major = order(abs(eigen.decomp[["values"]]), decreasing=TRUE)[1:k]
U = eigen.decomp[["vectors"]]

U.k = U[,major]
D.k = diag(eigen.decomp[["values"]][major])

TU = t(U) %*% T
TU.k = t(U.k) %*% T

QR.TU.k = qr(TU.k)
QR.TU = qr(TU)

Z.k = qr.Q(QR.TU.k)

X.mod = U.k %*% D.k %*% Z.k
XX = cbind(X.mod, T)

m1 = lm(income~XX-1)
delta = U.k %*% Z.k %*% as.matrix(m1$coef[1:3])
alpha = as.matrix(m1$coef[4:6])
