gwlars.adaptive.ssr <- function(bw, x, y, colocated, dist, s, verbose, prior.weights, gweight, target, mode=mode) {
    #Get the weights
    reps = length(colocated)
    loow = gweight(dist, bw)[-colocated]
    w <- prior.weights[-colocated] * loow
    
    #Build the model
    xx = as.matrix(x[-colocated,])
    yy = as.matrix(y[-colocated])

    m <- ncol(xx)
    n <- nrow(xx)
    one <- rep(1, n)
    meanx <- drop(one %*% xx)/n
    x.centered <- scale(xx, meanx, FALSE)         # first subtracts mean
    normx <- sqrt(drop(one %*% (x.centered^2)))
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
    
    lm.step = try(lm(yy~xs, weights=w))  # mle fit on standardized
    
    if(class(lm.step) == "try-error") { 
        cat(paste("Couldn't make a model for finding the SSR at bandwidth ", bw, "\n", sep=""))
        return(Inf)
    }

    beta.lm = lm.step$coeff[2:(m+1)]                    # mle except for intercept
    adapt.weight = abs(beta.lm)                         # weights for adaptive lasso
    for (k in 1:dim(x.centered)[2]) {
        if (!is.na(adapt.weight[k])) {
            xs[,k] = xs[,k] * adapt.weight[k]
        } else {
            xs[,k] = rep(0, dim(xs)[1])
            adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
        }
    }

    model = lars(x=xs, y=yy, weights=w, family=family, lambda=s)
    ll = model$lambda

    #Find lambda to minimize CV error
    xs.colocated = (x[colocated,] - meanx) * adapt.weight / normx
    predictions = predict(model, newx=matrix(xs.colocated, nrow=reps, ncol=dim(xs)[2]), s=ll, type='fit', mode=mode)[['fit']]
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
    s.optimal = ll[which.min(cv.error)]
    
    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xs, s=s.optimal, type='fit', mode=mode)[['fit']]    
    resid = sum(w * (yy - fitted)**2)
    return((resid-target)**2)
}
