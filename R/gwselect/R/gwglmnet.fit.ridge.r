gwglmnet.fit.ridge = function(x, y, coords, S, indx=NULL, loc, bw=NULL, dist=NULL, s=NULL, family, verbose, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, interact, precondition, tau=3) {
    if (!is.null(indx)) {
        colocated = which(round(coords[indx,1],5)==round(as.numeric(loc[1]),5) & round(coords[indx,2],5) == round(as.numeric(loc[2]),5))
    }
    else {
        colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))        
    }
    reps = length(colocated)

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } 
    gwr.weights = drop(gwr.weights)  

    if (!is.null(indx)) {
        gwr.weights = gwr.weights[indx]
    }

    #For interaction on location:
    if (interact) {
        newnames = vector()
        oldnames = colnames(x)
        for (l in 1:length(oldnames)) {
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[1], sep=""))
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[2], sep=""))
        }
        interacted = matrix(ncol=2*ncol(x), nrow=nrow(x))
        for (k in 1:ncol(x)) {
            interacted[,2*(k-1)+1] = x[,k]*(coords[,1]-loc[1,1])
            interacted[,2*k] = x[,k]*(coords[,2]-loc[1,2])
        }
        x.interacted = cbind(x, interacted)
        colnames(x.interacted) = c(oldnames, newnames)
    }



    xx = as.matrix(x)
    if (interact) {xx.interacted = as.matrix(x.interacted)}
    yy = as.matrix(y)
    w <- prior.weights * gwr.weights
    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
   
    xx = xx[weighted,]
    if (interact) {xx.interacted = xx.interacted[weighted,]}
    yy = as.matrix(yy[weighted])
    meany = mean(yy)
    yy = yy - meany
    w = w[weighted]
    colocated = which(gwr.weights[weighted]==1)

    if (interact) { 
        xxx = xx.interacted
    } else {
        xxx = xx
    }
    if (precondition==TRUE) {
        s = svd(xxx)
        F = s$u  %*% diag(1/sqrt(s$d**2 + tau))  %*%  t(s$u)
        xxx = F %*% xxx
        yyy = F %*% yyy
    }  

    one <- rep(1, nrow(xxx))
    meanx <- drop(one %*% xxx) / nrow(xxx)
    x.centered <- scale(xxx, meanx, FALSE)         # first subtracts mean
    normx <- sqrt(drop(one %*% (x.centered**2)))
    names(normx) <- NULL
    xs = x.centered
   
    for (k in 1:dim(x.centered)[2]) {
        if (normx[k]!=0) {
            xs[,k] = xs[,k] / normx[k]
        } else {
            xs[,k] = rep(0, dim(xs)[1])
            normx[k] = Inf #This should allow the lambda-finding step to work.
        }
    }

    n = length(yy)
    p = dim(xs)[2]

#    lambda.min = optimize(trace, x=xs, w=w, interval=c(0,1))$minimum
#print("lambda.min:")
#print(lambda.min)
#    lambda = seq(lambda.min*2, 0.1, length.out=10001)
    lambda = seq(0, 0.1, length.out=1001)
    ridge.pen = cv.glmnet(y=yy, x=xs, standardize=FALSE, intercept=TRUE, family=family, weights=w, alpha=0, lambda=lambda, nfolds=5)$lambda.1se
print("ridge.pen:")
print(ridge.pen)
    ridge = glmnet(y=yy, x=xs, standardize=FALSE, intercept=TRUE, family=family, weights=w, alpha=0, lambda=ridge.pen)

    coefs = as.matrix(coef(ridge))
    xtxswrI = solve(t(cbind(1,xs)) %*% diag(w) %*% cbind(1,xs) + diag(c(0,rep(ridge.pen, p))))
    tr = sum(w * diag(cbind(1,xs) %*% xtxswrI %*% t(cbind(1,xs))))
print("tr:")
print(tr)
    sig2 = sum(w*(yy - cbind(1,xs) %*% coef(ridge))**2) / (sum(w) - tr)
print("sig2:")
print(sig2)
    Vb = as.matrix(sig2 * xtxswrI)
print("coefs:")
print(coefs)
print("Vb:")
print(Vb)
    beta.resampled = t(mvrnorm(n=S, mu=coefs, Sigma=Vb))
    resampled = (cbind(1,xs) %*%  beta.resampled + meany)[colocated,]

    return(resampled)
}


trace = function(lambda, x, w) { 
    p = ncol(x)

    xtxswrI = solve(t(cbind(1,x)) %*% diag(w) %*% cbind(1,x) + diag(c(0,rep(lambda, p))))
    tr = sum(w * diag(cbind(1,x) %*% xtxswrI %*% t(cbind(1,x))))

    return((sum(w)-tr)**2)
}