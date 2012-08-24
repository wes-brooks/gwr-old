gwlars.ssr <- function(bw, x, y, colocated, dist, s, verbose, prior.weights, gweight, target, mode) {
    #Get the weights
    reps = length(colocated)
    loow = gweight(dist, bw)[-colocated]
    w <- prior.weights[-colocated] * loow
    
    #Build the model
    xx = as.matrix(x[-colocated,])
    yy = as.matrix(y[-colocated])

    m <- ncol(xx)
    n <- nrow(xx)
    
    model = lars(x=xx, y=yy, weights=w, family=family, lambda=s)
    ll = model$lambda

    #Find lambda to minimize CV error
    x.colocated = (x[colocated,] - meanx) * adapt.weight / normx
    predictions = predict(model, newx=matrix(x.colocated, nrow=reps, ncol=dim(xx)[2]), s=ll, type='fit', mode=mode)[['fit']]
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
    s.optimal = ll[which.min(cv.error)]
    
    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xx, s=s.optimal, type='fit', mode=mode)[['fit']]
    resid = sum(w * (yy - fitted)**2)
    return((resid-target)**2)
}

