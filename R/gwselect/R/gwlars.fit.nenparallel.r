gwlars.fit.nenparallel = function(x, y, coords, D, N=N, s, mode.select, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, mode, precondition=FALSE, AICc) {
    if (!is.null(fit.loc)) {
        coords.unique = fit.loc
    } else {
        coords.unique = unique(coords)
    }
    n = dim(coords.unique)[1]

    gwlars.object = list()
    models = list()

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(D, bw)    
    }       

    gweights = list()
    for (j in 1:nrow(gwr.weights)) {
        gweights[[j]] = as.vector(gwr.weights[j,])
    }

    if (verbose) {cat(paste('beta1:', beta1, ', beta2:', beta2, '\n', sep=''))}

    models = foreach(i=1:n, .packages=c('lars'), .errorhandling='remove') %dopar% {
        loc = coords.unique[i,]
        dist = D[i,]

        opt = optimize(gwlars.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, N=N, coords=coords, loc=loc, s=s,
            gweight=gweight, verbose=verbose, dist=dist, adapt=adapt, mode.select=mode.select,
            prior.weights=prior.weights, target=target, precondition=precondition)
        bandwidth = opt$minimum
		
		if (verbose) {
	        cat(paste("For i=", i, ", target: ", target, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", opt$objective, ".\n", sep=''))
        }
        return(gwlars.fit.inner(x=x, y=y, coords=coords, loc=loc, bw=bandwidth, dist=dist, N=N, s=s, mode.select=mode.select, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode=mode, precondition=precondition, AICc=AICc))
    }

    gwlars.object[['models']] = models
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['s.range']] = s

    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
