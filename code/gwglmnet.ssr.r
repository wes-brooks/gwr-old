gwglmnet.ssr <- function(bw, x, y, colocated, dist, s, verbose, family, prior.weights, gweight, type, target) {
    #Get the weights
    reps = length(colocated)
    loow = gweight(dist, bw)[-colocated]
    w <- prior.weights[-colocated] * loow
    
    #Build the model
    xx = as.matrix(x[-colocated,])
    yy = as.matrix(y[-colocated])

    if (family=='binomial') {model = glmnet(x=xx, y=cbind(1-yy,yy), weights=w, family=family, lambda=s)}
    else {model = glmnet(x=xx, y=yy, weights=w, family=family, lambda=s)}
    ll = model$lambda

    #Find lambda to minimize CV error
    predictions = predict(model, newx=matrix(x[colocated,], nrow=reps, ncol=dim(xx)[2]), s=ll, type='response', )
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
    s.optimal = ll[which.min(cv.error)]
    
    #Get the residuals at this choice of s (Poisson-specific for now):
    fitted = predict(model, newx=xx, s=s.optimal, type='response')    
    if (family=='poisson') pearson.resid = sum(w * (yy - fitted)**2/fitted)
    if (family=='binomial') pearson.resid = sum(w * (yy - fitted)**2 / (fitted*(1-fitted)))

    (abs(pearson.resid-target))**2
}