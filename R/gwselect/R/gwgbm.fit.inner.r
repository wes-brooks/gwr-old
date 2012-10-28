gwgbm.fit.inner = function(x, y, family, coords, loc, bw=NULL, dist=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, precondition=FALSE) {
    colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)
    }
    
    w <- prior.weights * gwr.weights   
    loow = w[-colocated] 
    reps = length(colocated)        

    if (sum(gwr.weights[-colocated])==0) { return(list(loss.local=Inf, resid=Inf)) }   
       
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

    
	meanx = rep(0, ncol(x))
	adapt.weight = rep(1, ncol(x))
	normx = rep(1, ncol(x))

	xs=xx
	predx = matrix(x[colocated,], nrow=reps, ncol=ncol(x)) 


    xfit = data.matrix(xs)
    yfit = data.matrix(yy)
    fitdata = data.frame(xfit, yfit)

    if (family=='binomial' && (abs(sum(yfit*loow)-sum(w))<1e-4 || sum(yfit*loow)<1e-4)) {
        cat(paste("Abort. i=", i, ", weighted sum=", sum(yfit*loow), ", sum of weights=", sum(loow), "\n", sep=''))
        return(list(model=NULL, cv.error=0, s=Inf, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=Inf))
    } else if (family=='binomial') {
        model = gbm(yfity~., data=fitdata, distribution=family, weights=w, n.trees=5000, bag.fraction=0.5, interaction.depth=5, n.minobsinnode=10)
        predictions = predict(model, newx=predx, type='response')
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(model[['lambda']])), nrow=reps, ncol=length(model[['lambda']]))))
        s.optimal = model[['lambda']][which.min(cv.error)]
        V = function(mu) {mu*(1-mu)}
    } else if (family=='poisson') {
        model = glmnet(x=xfit, y=yfit, family=family, weights=loow, lambda=s)
        predictions = predict(model, newx=predx, type='response')
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(model[['lambda']])), nrow=reps, ncol=length(model[['lambda']]))))
        s.optimal = model[['lambda']][which.min(cv.error)]
        V = function(mu) {mu}
    }

    #Get the coefficients:
    coef = predict(model, type='coefficients', s=s.optimal, mode='lambda')

    coef[2:length(coef)] = coef[2:length(coef)] * adapt.weight/normx
    coef[1] = coef[1] - sum(coef[2:length(coef)] * meanx)

    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xfit, s=s.optimal, type='response')
    resid = sqrt(w[-colocated]) * (yfit - fitted) / sqrt(V(fitted))    

    return(list(model=model, cv.error=cv.error, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=resid, coef=coef))
}