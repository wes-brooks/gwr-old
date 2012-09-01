gwglmnet.fit.inner = function(x, y, family, coords, loc, bw=NULL, dist=NULL, s=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE) {
    colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)
    }
    
    w <- prior.weights * gwr.weights   
    loow = w[-colocated] 
    reps = length(colocated)        

    if (sum(gwr.weights[-colocated])==0) { return(list(cv.error=Inf, resid=Inf)) }   
       
    xx = as.matrix(x[-colocated,])
    yy = as.matrix(y[-colocated])
    
    if (precondition==TRUE) {
        s = svd(xx)
        F = s$u  %*% diag(1/s$d)  %*%  t(s$u)
        xx = F %*% xx
        yy = F %*% yy
    }

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
        
        glm.step = try(glm(yy~xs, weights=loow, family=family))  # mle fit on standardized
    
        if(class(glm.step) == "try-error") { 
            cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
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
        predx = (x[colocated,] - meanx) * adapt.weight / normx
    } else {
        meanx = rep(0, ncol(x))
        adapt.weight = rep(1, ncol(x))
        normx = rep(1, ncol(x))

        xs=xx
        predx = x[colocated,]
    }

    xfit = xs
    yfit = yy

    if (family=='binomial' && (abs(sum(yfit*loow)-sum(w))<1e-4 || sum(yfit*loow)<1e-4)) {
        cat(paste("Abort. i=", i, ", weighted sum=", sum(yfit*loow), ", sum of weights=", sum(loow), "\n", sep=''))
        return(list(model=NULL, cv.error=0, s=Inf, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=Inf))
    } else if (family=='binomial') {
        model = glmnet(x=xfit, y=cbind(1-yfit, yfit), family=family, weights=loow, lambda=s)
        predictions = predict(model, newx=predx, type='response')
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(model[['lambda']])), nrow=reps, ncol=length(model[['lambda']]))))
        s.optimal = model[['lambda']][which.min(cv.error)]
        V = function(mu) {mu*(1-mu)}
    } else {
        model = glmnet(x=xfit, y=yfit, family=family, weights=loow, lambda=s)
        predictions = predict(model, newx=predx, type='response')
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(model[['lambda']])), nrow=reps, ncol=length(model[['lambda']]))))
        s.optimal = model[['lambda']][which.min(cv.error)]
        V = function(mu) {mu}
    }

    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xfit, s=s.optimal, type='response')
    resid = sqrt(w[colocated]) * (yfit - fitted) / sqrt(V(fitted))    

    return(list(model=model, cv.error=cv.error, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=resid))
}