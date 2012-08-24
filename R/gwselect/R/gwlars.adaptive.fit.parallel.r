gwlars.nen.fit.parallel = function(x, y, coords, D, s, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE) {
    #Fit the GWLARS model (non-adaptive algorithm)
    coords.unique = unique(coords)

    gwlars.object = foreach(i=1:n, .packages='lars', .errorhandling='remove') %dopar% {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
        dist = D[i,]

        opt = optimize(gwlars.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, colocated=colocated, s=s,
            gweight=gweight, verbose=verbose, dist=dist,
            prior.weights=prior.weights, target=target, type=type)
        bandwidth = opt$minimum

        cat(paste("For i=", i, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", sqrt(opt$objective), ".\n", sep=''))

        loow = gweight(D[i,-colocated], bandwidth)        
        prior.loow = prior.weights[-colocated]
        w <- prior.loow * loow
        reps = length(colocated)      


        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
    
        
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
       

        xs.colocated = (x[colocated,] - meanx) * adapt.weight / normx
        model = glmnet(x=xs, y=cbind(1-yy, yy), weights=w, family=family, lambda=s)
        ll = model$lambda
        predictions = predict(model, newx=matrix(xs.colocated, nrow=reps, ncol=dim(xx)[2]), s=ll, type='fit')[['fit']]
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
        s.optimal = ll[which.min(cv.error)]
        
        if (verbose) { cat(paste(i, "\n", sep='')) }         
        list(model=model, cv.error=cv.error, s=s.optimal, index=i)
    }

    print("returning from gwlars.nen.adaptive.fit.parallel")
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}



#        colocated = which(coords$x==coords.unique$x[i] & coords$y==coords.unique$y[i])
#        loow = weight.matrix[i,-colocated]
#        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
#        reps = length(colocated)        
#        w.sqrt <- diag(rep(sqrt(loow), reps))
        
#        model[[i]] = lars(x=w.sqrt %*% as.matrix(x[-colocated,]), y=w.sqrt %*% as.matrix(y[-colocated]))
#        predictions = predict(model[[i]], newx=matrix(x[colocated,], nrow=reps, ncol=dim(x)[2]), s=s, type='fit', mode=mode)[['fit']]
#        cv.error[[i]] = colSums(abs(matrix(predictions - as.matrix(y[colocated]), nrow=reps, ncol=length(s))))
#        s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])