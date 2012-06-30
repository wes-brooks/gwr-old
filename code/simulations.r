size = c(20, 30, 40, 50)

for (N in size) {
    coord = seq(0, 1, length.out=N)
    grid = matrix(rnorm(N**2, mean=rep(ifelse(coord<=0.5, 0, 2), N), sd=0.5), N, N)
    rownames(grid) = seq(0, 1, length.out=N)
    colnames(grid) = seq(0, 1, length.out=N)
}

#
X1 = matrix(rnorm(N**2, mean=0, sd=1), N, N)
B1 = matrix(rep(ifelse(coord<=0.5, 0, 2), N), N, N)
Y = X1 * B1 + rnorm(N**2, 0, 0.1)
Z = rnorm(N**2, 0, 1)

#
loc.x = rep(seq(0, 1, length.out=N), each=N)
loc.y = rep(seq(0, 1, length.out=N), times=N)
sim = data.frame(Y=as.vector(Y),X=as.vector(X1), Z, loc.x, loc.y)

#
n = dim(sim)[1]
Xmat = matrix(rep(loc.x,n), n,n)
Ymat = matrix(rep(loc.y,n), n,n)
D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)

