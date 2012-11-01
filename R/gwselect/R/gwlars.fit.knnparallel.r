gwlars.fit.knnparallel = function(x, y, coords, indx, fit.loc, D, N=N, s, mode.select, tuning, predict, shrink, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, mode, precondition=FALSE) {
    if (!is.null(fit.loc)) {
        coords.unique = unique(fit.loc)
        n.loc = dim(coords.unique)[1]
    } else {
        coords.unique = unique(coords)
        n.loc = dim(coords.unique)[1]
    }

    n = dim(coords.unique)[1]
    gwlars.object = list()

    prior.weights = drop(prior.weights)
    max.weights = rep(1, length(prior.weights))
    total.weight = sum(max.weights * prior.weights)

    models = foreach(i=1:n, .packages=c('lars'), .errorhandling='remove') %dopar% {
        loc = coords.unique[i,]
        dist = drop(D[i,])

        opt = optimize(gwlars.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, coords=coords, indx=indx, loc=loc, s=s,
            gweight=gweight, verbose=verbose, dist=dist, adapt=adapt, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ".\n", sep=''))
        return(gwlars.fit.inner(x=x, y=y, coords=coords, indx=indx, loc=loc, bw=bandwidth, dist=dist, N=N, s=s, mode.select=mode.select, tuning=tuning, predict=predict, shrink=shrink, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode=mode, precondition=precondition))
    }

    gwlars.object[['models']] = models
    gwlars.object[['mode']] = mode

    if (tuning) {
    } else if (predict) {
    } else {
        gwlars.object[['coords']] = coords
        gwlars.object[['s.range']] = s
    }

    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
