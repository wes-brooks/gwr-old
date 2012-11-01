gwglmnet.fit.knnparallel = function(x, y, family, coords, fit.loc, indx, D, s, verbose, prior.weights, tuning, predict, gweight, mode, mode.select, shrink, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE, N) {
    if (!is.null(fit.loc)) {
        coords.unique = unique(fit.loc)
        n.loc = dim(coords.unique)[1]
    } else {
        coords.unique = unique(coords)
        n.loc = dim(coords.unique)[1]
    }

    n = dim(prior.weights)[1]
    gwglmnet.object = list()

    prior.weights = drop(prior.weights)
    max.weights = rep(1, n)
    total.weight = sum(max.weights * prior.weights)

    models = foreach(i=1:n.loc, .packages=c('glmnet'), .errorhandling='remove') %dopar% {
        loc = coords.unique[i,]
        dist = drop(D[i,])

        opt = optimize(gwglmnet.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, coords=coords, loc=loc, indx=indx,
            gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum
        
        m = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, tuning=tuning, predict=predict, indx=indx, loc=loc, bw=bandwidth, dist=dist, s=s, verbose=verbose, mode=mode, mode.select=mode.select, shrink=shrink, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, N=N)
        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ", loss=", m[['loss.local']], ".\n", sep=''))
        return(m)
    }

    gwglmnet.object[['models']] = models

    if (tuning) {
    } else if (predict) {
    } else {
        gwglmnet.object[['coords']] = coords
        gwglmnet.object[['s.range']] = s
    }

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
