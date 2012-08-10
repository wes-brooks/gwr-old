gwglmnet.nen.fit.parallel = function(x, y, coords, D, s, verbose, family, prior.weights, gweight, target, beta1, beta2, type='pearson', tol=1e-25, longlat=FALSE) {
    #Fit the gwglmnet model (non-adaptive algorithm)
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]


    gwglmnet.object = foreach(i=1:n, .packages='glmnet', .errorhandling='remove') %dopar% {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
        dist = D[i,]

        opt = optimize(gwglmnet.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, colocated=colocated, s=s,
            gweight=gweight, verbose=verbose, dist=dist,
            prior.weights=prior.weights, family=family, target=target, type=type)
        bandwidth = opt$minimum
        
        cat(paste("For i=", i, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", sqrt(opt$objective), ".\n", sep=''))

        loow = gweight(D[i,-colocated], bandwidth)        
        prior.loow = prior.weights[-colocated]
        w <- prior.loow * loow
        reps = length(colocated)        

        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
            
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
        

        if (family=='binomial' && (abs(sum(yy*w)-sum(w))<1e-4 || sum(yy*w)<1e-4)) {
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yy*w), ", sum of weights=", sum(w), "\n", sep=''))
            model = NULL
            cv.error = 0
            s.optimal = max(s)
        } else if (family=='binomial') {
            model = glmnet(x=xx, y=cbind(1-yy, yy), weights=w, family=family, lambda=s)
            predictions = predict(model, newx=matrix(x[colocated,], nrow=reps, ncol=dim(xx)[2]), s=s, type='response')
            cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(s)), nrow=reps, ncol=length(s))))
            s.optimal = s[which.min(cv.error)]
        } else {
            model = glmnet(x=xx, y=yy, weights=w, family=family, lambda=s)
            ll = model$lambda
            predictions = predict(model, newx=matrix(x[colocated,], nrow=reps, ncol=dim(xx)[2]), s=ll, type='response')
            cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
            s.optimal = ll[which.min(cv.error)]
        }
        
        if (verbose) { cat(paste(i, "\n", sep='')) }         

        list(model=model, cv.error=cv.error, s=s.optimal, index=i)
    }

    print("returning from gwglmnet.nen.fit.parallel")

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
