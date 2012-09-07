gwglmnet.fit.nenparallel = function(x, y, family, coords, D, s, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwglmnet.object = list()

    if (verbose) {cat(paste('beta1:', beta1, ', beta2:', beta2, '\n', sep=''))}

    models = foreach(i=1:n, .packages=c('glmnet'), .errorhandling='remove') %dopar% {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwglmnet.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, coords=coords, loc=loc, s=s,
            gweight=gweight, verbose=verbose, dist=dist, adapt=adapt, family=family,
            prior.weights=prior.weights, target=target, precondition=precondition)
        bandwidth = opt$minimum

        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ".\n", sep=''))
        return(gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bandwidth, dist=dist, s=s, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition))
    }

    gwglmnet.object[['models']] = models
    gwglmnet.object[['coords']] = coords
    gwglmnet.object[['s.range']] = s

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
