gwlars.fit.fixedbwparallel = function(x, y, coords, indx, fit.loc, bw, D=NULL, N, s=NULL, mode.select, tuning, predict, simulation, oracle, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE, interact, shrunk.fit, AICc) {
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

    models = foreach(i=1:n, .packages=c('lars'), .errorhandling='remove') %dopar% {
        #Fit one location's model here
        loc = coords.unique[i,]
        gw = gweights[[i]]

        if (is.null(oracle)) {
            m = gwlars.fit.inner(x=x, y=y, bw=bw, coords=coords, loc=loc, indx=indx, N=N, s=s, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, interact=interact, shrunk.fit=shrunk.fit, AICc=AICc)
        } else {
            m = gwlars.fit.oracle(x=x, y=y, family='gaussian', bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, interact=interact, AICc=AICc)
        }
        
        if (verbose) {
        	cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bw,3), "; loss=", round(m[['loss.local']],3), "; s=", m[['s']], "; sigma2=", round(m[['sigma2']],3), "; nonzero=", paste(m[['nonzero']], collapse=", "), ".\n", sep=''))
        }
        return(m)
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
