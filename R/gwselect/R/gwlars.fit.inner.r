gwlars.fit.inner = function(x, y, coords, loc, bw=NULL, dist=NULL, s=NULL, mode.select='', verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, mode, precondition=FALSE) {
    colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))
    reps = length(colocated)

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } else {
        gwr.weights = gwr.weights
    }      

    if (sum(gwr.weights)==0) { return(list(cv.error=Inf, resid=Inf)) }   
      
    if (mode.select=='CV') { 
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
        w <- prior.weights[-colocated] * gwr.weights[-colocated]   
    } else {
        xx = as.matrix(x)
        yy = as.matrix(y)
        w <- prior.weights * gwr.weights
    }
    
    sqrt.w <- diag(sqrt(w))       
    
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
        
        lm.step = try(lm(yy~xs, weights=w))  # mle fit on standardized
    
        if(class(lm.step) == "try-error") { 
            cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
            return(Inf)
        }
        
        beta.lm = lm.step$coeff[2:(m+1)]                   # mle except for intercept
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
        adapt.weight = rep(1, ncol(x))
        normx = rep(1, ncol(x))

        xs=xx
        predx = x[colocated,]
    }

    xfit = sqrt.w %*% xs
    yfit = sqrt.w %*% yy
    model = lars(x=xfit, y=yfit, type='lar', normalize=FALSE, intercept=TRUE)
    #ll = model$lambda
    nsteps = length(model$lambda) + 1

    if (mode.select=='CV') {
        predx = matrix(predx, reps, dim(xs)[2])
        #predictions = predict(model, newx=predx, s=ll, type='fit', mode='lambda')[['fit']]
        predictions = predict(model, newx=predx, type='fit')[['fit']]
        #cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
        loss = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=nsteps), nrow=reps, ncol=nsteps)))
        #s.optimal = ll[which.min(cv.error)]
    } else if (mode.select=='AIC') {
        coef = predict(model, type='coefficients')[['coefficients']]
        df = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {sum(abs(x)>0)}) + 1
        fitted = matrix(predict(model, newx=xfit, type='fit', mode='lambda')[['fit']], n, nsteps)
        s2 = (fitted[,nsteps] - yfit)**2 / sum(w)
        loss = as.vector(apply(fitted, 2, function(z) {sum((z - yfit)**2)})/(s2*sum(w)) + 2*df/sum(w))
        #s.optimal = which.min(    
    } else if (mode.select=='BIC') {
        coef = predict(model, type='coefficients')[['coefficients']]
        df = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {sum(abs(x)>0)}) + 1
        fitted = matrix(predict(model, newx=xfit, type='fit', mode='step')[['fit']], n, nsteps)
        s2 = (fitted[,nsteps] - yfit)**2 / sum(w)
        loss = as.vector(apply(fitted, 2, function(z) {sum((z - yfit)**2)})/(s2*sum(w)) + log(sum(w))*df/sum(w)) 
    }

    print(loss)
    s.optimal = numeric(which.min(loss))
    print(s.optimal)
    
    #Get the coefficients:
    coef = predict(model, type='coefficients', s=s.optimal, mode='step')[['coefficients']]
    coef = Matrix(coef, ncol=1)
    rownames(coef) = colnames(x)

    intercept = predict(model, type='fit', s=s.optimal, mode='step', newx=matrix(0,1,nrow(coef)))[['fit']]

    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xfit, s=s.optimal, type='fit', mode='step')[['fit']]
    resid = yfit - fitted
    
    return(list(model=model, cv.error=loss, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=resid, coef=coef, intercept=intercept))
}