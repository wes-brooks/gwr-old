size = c(20, 30, 40, 50)

for (N in size) {
    coord = seq(0, 1, length.out=N)
    grid = matrix(rnorm(N**2, mean=rep(ifelse(coord<=0.5, 0, 2), N), sd=0.5), N, N)
    rownames(grid) = seq(0, 1, length.out=N)
    colnames(grid) = seq(0, 1, length.out=N)
}

#population for weights
pop = rpois(N**2, 400)

#
X1 = matrix(rnorm(N**2, mean=0, sd=1), N, N)
B1 = matrix(rep(ifelse(coord<=0.5, 0, 2), N), N, N)

X2 = matrix(rnorm(N**2, mean=0, sd=1), N, N)
B2 = matrix(rep(1-coord, N), N, N)

eta = X1*B1 + X2*B2 + rnorm(N**2, 0, 0.1)
Z = rnorm(N**2, 0, 1)
p = exp(eta) / (1+exp(eta))
Y = rbinom(N**2, pop, p) / pop

#
loc.x = rep(seq(0, 1, length.out=N), each=N)
loc.y = rep(seq(0, 1, length.out=N), times=N)
sim = data.frame(Y=as.vector(Y), X1=as.vector(X1), X2=as.vector(X2), Z, loc.x, loc.y)

#
n = dim(sim)[1]
Xmat = matrix(rep(loc.x,n), n,n)
Ymat = matrix(rep(loc.y,n), n,n)
D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)

