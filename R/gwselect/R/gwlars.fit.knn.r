gwlars.fit.knn = function(x, y, coords, indx, fit.loc, D, N=N, s, mode.select, tuning, predict, simulation, verbose, interact, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE, oracle, shrunk.fit, AICc) {
    if (!is.null(fit.loc)) {
        coords.unique = fit.loc
    } else {
        coords.unique = unique(coords)
    }
    n = dim(coords.unique)[1]

    gwlars.object = list()
    models = list()

    prior.weights = drop(prior.weights)
    max.weights = rep(1, length(prior.weights))
    total.weight = sum(max.weights * prior.weights)

    for (i in 1:n) {
        loc = coords.unique[i,]
        dist = drop(D[i,])

        opt = optimize(gwlars.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, coords=coords, indx=indx, loc=loc, s=s,
            gweight=gweight, verbose=verbose, dist=dist, adapt=adapt, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum

        if (is.null(oracle)) {
            models[[i]] = gwlars.fit.inner(x=x, y=y, bw=bandwidth, coords=coords, loc=loc, indx=indx, N=N, s=s, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, dist=dist, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, interact=interact, shrunk.fit=shrunk.fit, AICc=AICc)
        } else {
            models[[i]] = gwselect.fit.oracle(x=x, y=y, bw=bandwidth, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, dist=dist, prior.weights=prior.weights, gweight=gweight, interact=interact, AICc=AICc)
        }
        
        if (verbose) {
	        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); target=", target, "; bw=", bandwidth, "; tolerance=", target/1000, "; loss=", models[[i]][['loss.local']], ".\n", sep=''))
		}
    }

    gwlars.object[['models']] = models

    if (tuning) {
    } else if (predict) {
    } else if (simulation) {
    } else {
        gwlars.object[['coords']] = coords
        gwlars.object[['s.range']] = s
    }

    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
