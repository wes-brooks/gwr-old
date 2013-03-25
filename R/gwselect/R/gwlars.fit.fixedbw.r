gwlars.fit.fixedbw = function(x, y, coords, indx, fit.loc, bw, D=NULL, N, s=NULL, mode.select, tuning, predict, simulation, oracle, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE, interact) {
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

    for (i in 1:n) {
        #Fit one location's model here
        loc = coords.unique[i,]
        gw = gweights[[i]]

        if (is.null(oracle)) {
            models[[i]] = gwlars.fit.inner(x=x, y=y, bw=bw, coords=coords, loc=loc, indx=indx, N=N, s=s, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, interact=interact)
        } else {
            models[[i]] = gwlars.fit.oracle(x=x, y=y, bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, interact=interact)
        }

        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bw,3), "; loss=", round(models[[i]][['loss.local']],12), "; s=", models[[i]][['s']], "; sigma2=", round(models[[i]][['sigma2']],3), "; nonzero=", paste(models[[i]][['nonzero']], collapse=","), "; weightsum=", round(models[[i]][['weightsum']],3), ".\n", sep=''))
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
