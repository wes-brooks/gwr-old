gwglmnet.fit.knn = function(x, y, family, coords, fit.loc, D, s, verbose, prior.weights, tuning, predict, indx, gweight, mode, mode.select, shrink, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE, N) {
    
    n = dim(prior.weights)[1]

    if (!is.null(fit.loc)) {
        coords.unique = unique(fit.loc)
        n.loc = dim(coords.unique)[1]
    } else {
        coords.unique = unique(coords)
        n.loc = dim(coords.unique)[1]
    }
    gwglmnet.object = list()
    models = list()

    prior.weights = drop(prior.weights)
    max.weights = rep(1, n)
    total.weight = sum(max.weights * prior.weights)

    for (i in 1:n.loc) {
        loc = coords.unique[i,]
        dist = drop(D[i,])

        opt = optimize(gwglmnet.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, coords=coords, loc=loc, indx=indx,
            gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, tuning=tuning, predict=predict, indx=indx, bw=bandwidth, dist=dist, s=s, verbose=verbose, mode=mode, mode.select=mode.select, shrink=shrink, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, N=N)
        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ", loss=", models[[i]][['loss.local']], ".\n", sep=''))
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
