gwlars.fit.nen = function(x, y, coords, D, s, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, mode) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwlars.object = list()
    models = list()

    if (verbose) {cat(paste('beta1:', beta1, ', beta2:', beta2, '\n', sep=''))}

    for (i in 1:n) {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwlars.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, coords=coords, loc=loc, s=s,
            gweight=gweight, verbose=verbose, dist=dist, adapt=adapt,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

        models[[i]] = gwlars.fit.inner(x=x, y=y, coords=coords, loc=loc, bw=bandwidth, dist=dist, s=s, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode=mode)

        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ".\n", sep=''))
    }

    gwlars.object[['models']] = models
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['s.range']] = s

    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
