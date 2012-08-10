gwglmnet.adaptive.ssr <- function(bw, x, y, colocated, dist, s, verbose, family, prior.weights, gweight, type, target) {
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
    
    glm.step = try(glm(yy~xs, family=family, weights=w))  # mle fit on standardized
    
    if(class(glm.step) == "try-error") { 
        cat(paste("Couldn't make a model for finding the SSR at bandwidth ", bw, "\n", sep=""))
        #glm.step[[i]] = out.glm = glm.step[[i-1]]
        return(Inf)
    }

    beta.glm = glm.step$coeff[2:(m+1)]                    # mle except for intercept
    adapt.weight = abs(beta.glm)                        # weights for adaptive lasso
    for (k in 1:dim(x.centered)[2]) {
        if (!is.na(adapt.weight[k])) {
            xs[,k] = xs[,k] * adapt.weight[k]
        } else {
            xs[,k] = rep(0, dim(xs)[1])
            adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
        }
    }

    if (family=='binomial') {model = glmnet(x=xs, y=cbind(1-yy,yy), weights=w, family=family, lambda=s)}
    else {model = glmnet(x=xs, y=yy, weights=w, family=family, lambda=s)}
    ll = model$lambda

    #Find lambda to minimize CV error
    xs.colocated = (x[colocated,] - meanx) * adapt.weight / normx
    predictions = predict(model, newx=matrix(xs.colocated, nrow=reps, ncol=dim(xs)[2]), s=ll, type='response', )
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
    s.optimal = ll[which.min(cv.error)]
    
    #Get the residuals at this choice of s (Poisson-specific for now):
    fitted = predict(model, newx=xs, s=s.optimal, type='response')    
    if (family=='poisson') pearson.resid = sum(w * (yy - fitted)**2/fitted)
    if (family=='binomial') pearson.resid = sum(w * (yy - fitted)**2 / (fitted*(1-fitted)))

    #cat(paste("bw: ", bw, ", cv error: ", ((abs(pearson.resid-target))**2), "\n", sep=""))
    (abs(pearson.resid-target))**2
}