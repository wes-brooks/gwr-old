gwglmnet.fit.fixedbw = function(x, y, family, coords, fit.loc=NULL, oracle, bw, D=NULL, s=NULL, verbose=FALSE, mode.select, gwr.weights=NULL, indx=NULL, prior.weights=NULL, tuning=FALSE, predict=FALSE, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE, alpha, simulation, interact, N, shrunk.fit, AICc) {
    if (!is.null(fit.loc)) {
        coords.unique = fit.loc
    } else {
        coords.unique = unique(coords)
    }
    n = dim(coords.unique)[1]

    gwglmnet.object = list()
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
	        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, s=s, verbose=verbose, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, predict=predict, tuning=tuning, simulation=simulation, alpha=alpha, interact=interact, N=N, shrunk.fit=shrunk.fit, AICc=AICc)
		} else {
            models[[i]] = gwselect.fit.oracle(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, AICc=AICc)
        }
        
        if (verbose) {
	        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bw, 3), "; loss=", round(tail(models[[i]][['loss.local']],1),3), "; s=", models[[i]][['s']], "; sigma2=", round(tail(models[[i]][['sigma2']],1),3), "; nonzero=", paste(models[[i]][['nonzero']], collapse=","), "; weightsum=", round(models[[i]][['weightsum']],3), ".\n", sep=''))
		}
    }

    gwglmnet.object[['models']] = models
	
	if (tuning) {
    } else if (predict) {
    } else if (simulation) {
    } else {
	    gwglmnet.object[['coords']] = coords
    	gwglmnet.object[['s.range']] = s
	}

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
