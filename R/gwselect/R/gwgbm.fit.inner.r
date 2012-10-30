gwgbm.fit.inner = function(x, y, family, coords, loc, bw=NULL, dist=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, precondition=FALSE) {
    colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))
    reps = length(colocated)
    
    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)
    }
    
    w <- prior.weights * gwr.weights         

    if (sum(gwr.weights[-colocated])==0) { return(list(loss.local=Inf, resid=Inf)) }   
       
    xx = as.matrix(x)
    yy = as.matrix(y)
    
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


	predx = matrix(x[colocated,], nrow=reps, ncol=ncol(x)) 


    xfit = data.matrix(xx)
    yfit = data.matrix(yy)

    #if (family=='binomial' && (abs(sum(yfit*loow)-sum(w))<1e-4 || sum(yfit*loow)<1e-4)) {
    #    cat(paste("Abort. i=", i, ", weighted sum=", sum(yfit*loow), ", sum of weights=", sum(loow), "\n", sep=''))
    #    return(list(model=NULL, cv.error=0, s=Inf, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=Inf))
    #} else 
    if (family=='binomial') {    
	    p = c(yfit, (1-yfit))
	    ww = rep(w, 2) * p
    	yfit2 = as.matrix(c(rep(1, dim(xfit)[1]), rep(0, dim(xfit)[1])))
    	xfit2 = rbind(xfit, xfit)
	    fitdata = data.frame(xfit2, yfit2)
	        
        model = gbm(yfit2~., data=fitdata, distribution='bernoulli', weights=ww, n.trees=10000, bag.fraction=0.5, interaction.depth=5, n.minobsinnode=10)
        s.optimal = gbm.perf(model, method="OOB", plot.it=FALSE)
        fitted = predict(model, newdata=as.data.frame(xfit), n.trees=s.optimal)
		loss.local = sum((w*(fitted - yy)**2 / (fitted*(1-fitted)))[colocated]) 
    } else if (family=='poisson') {
        model = glmnet(x=xfit, y=yfit, family=family, weights=loow, lambda=s)
        predictions = predict(model, newx=predx, type='response')
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(model[['lambda']])), nrow=reps, ncol=length(model[['lambda']]))))
        s.optimal = model[['lambda']][which.min(cv.error)]
        V = function(mu) {mu}
    }

    #Get the coefficients:
    coef = Matrix(rep(NA, ncol(x)))

    return(list(model=model, loss.local=loss.local, s=s.optimal, loc=loc, bw=bw, resid=resid, coef=coef))
}