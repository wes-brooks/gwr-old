gwgbm.fit.knnparallel = function(x, y, family, coords, D, s,mode.select=mode.select, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwgbm.object = list()

    max.weights = rep(1, n)
    total.weight = sum(max.weights * prior.weights)

    models = foreach(i=1:n, .packages=c('glmnet'), .errorhandling='remove') %dopar% {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwgbm.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, coords=coords, loc=loc, mode.select=mode.select,
            gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ".\n", sep=''))
        return(gwgbm.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bandwidth, dist=dist, s=s, mode.select=mode.select, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition))
    }

    gwgbm.object[['models']] = models
    gwgbm.object[['coords']] = coords
    gwgbm.object[['s.range']] = s

    class(gwgbm.object) = 'gwglmnet.object'
    return(gwgbm.object)
}
