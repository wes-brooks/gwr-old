gwlars.fit.inner = function(x, y, coords, loc, bw=NULL, dist=NULL, s=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, mode) {
    colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist[-colocated], bw)     
    } else {
        gwr.weights = gwr.weights[-colocated]
    }
    
    prior.loow = prior.weights[-colocated]
    w <- prior.loow * gwr.weights        
    reps = length(colocated)
    sqrt.w <- diag(sqrt(w))             

    if (sum(gwr.weights)==0) { return(list(cv.error=Inf, resid=Inf)) }   
       
    xx = as.matrix(x[-colocated,])
    yy = as.matrix(y[-colocated])
    
    m <- ncol(xx)
    n <- nrow(xx)

    if (adapt==TRUE) {
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
            cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
            return(Inf)
        }
    
        beta.lm = lm.step$coeff[2:(m+1)]                    # mle except for intercept
        adapt.weight = abs(beta.lm)                        # weights for adaptive lasso
        for (k in 1:dim(x.centered)[2]) {
            if (!is.na(adapt.weight[k])) {
                xs[,k] = xs[,k] * adapt.weight[k]
            } else {
                xs[,k] = rep(0, dim(xs)[1])
                adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
            }
        }
        
        predx = (x[colocated,] - meanx) * adapt.weight / normx
        
    } else {
        meanx = rep(0, ncol(x))
        adapt.weight[k] = rep(1, ncol(x))
        normx[k] = rep(1, ncol(x))

        xs=xx
        predx = x[colocated,]
    }

    xfit = sqrt.w %*% xs
    yfit = sqrt.w %*% yy
    model = lars(x=xfit, y=yfit, type='lar', normalize=FALSE)
    ll = model$lambda

    predictions = predict(model, newx=predx, s=ll, type='fit', mode='lambda')[['fit']]
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
    s.optimal = ll[which.min(cv.error)]

    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xfit, s=s.optimal, type='fit', mode='lambda')[['fit']]
    resid = yfit - fitted
    
    if (verbose) { cat(paste(i, "\n", sep='')) }         
    return(list(model=model, cv.error=cv.error, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=resid))
}