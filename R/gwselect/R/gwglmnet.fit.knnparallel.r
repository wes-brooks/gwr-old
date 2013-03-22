gwglmnet.fit.knnparallel = function(x, y, family, coords, fit.loc, indx, D, s, verbose, mode.select, prior.weights, tuning, predict, simulation, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt, precondition=FALSE, N, interact, alpha) {
    if (!is.null(fit.loc)) {
        coords.unique = unique(fit.loc)
    } else {
        coords.unique = unique(coords)
    }
    n = dim(coords.unique)[1]

    gwglmnet.object = list()

    prior.weights = drop(prior.weights)
    max.weights = rep(1, length(prior.weights))
    total.weight = sum(max.weights * prior.weights)

    models = foreach(i=1:n, .packages=c('glmnet'), .errorhandling='remove') %dopar% {
        loc = coords.unique[i,]
        dist = drop(D[i,])

        opt = optimize(gwglmnet.knn, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, coords=coords, loc=loc, indx=indx,
            gweight=gweight, verbose=verbose, dist=dist, total.weight=total.weight,
            prior.weights=prior.weights, target=target)
        bandwidth = opt$minimum
        
        if (is.null(oracle)) {
	        m = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, indx=indx, loc=loc, bw=bandwidth, dist=dist, s=s, verbose=verbose, mode.select=mode.select, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, N=N, interact=interact, alpha=alpha)
        } else {
            m = gwlars.fit.oracle(x=x, y=y, bw=bandwidth, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, dist=dist, prior.weights=prior.weights, gweight=gweight, interact=interact)
        }
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
