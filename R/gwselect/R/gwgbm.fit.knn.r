gwgbm.fit.knn = function(x, y, family, coords, D, s, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwgbm.object = list()
    models = list()

    max.weights = rep(1, n)
    total.weight = sum(max.weights * prior.weights)

    for (i in 1:n) {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwgbm.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, coords=coords, loc=loc,
            gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

        models[[i]] = gwgbm.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bandwidth, dist=dist, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, precondition=precondition)
        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ".\n", sep=''))
    }

    gwgbm.object[['models']] = models
    gwgbm.object[['coords']] = coords
    gwgbm.object[['s.range']] = s

    class(gwgbm.object) = 'gwgbm.object'
    return(gwgbm.object)
}
