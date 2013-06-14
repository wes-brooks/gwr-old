gwglmnet.fit.fixedbwparallel = function(x, y, family, coords, fit.loc=NULL, indx=NULL, bw, D=NULL, s=NULL, oracle, verbose=FALSE, mode.select, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat, adapt, precondition, alpha, simulation, tuning, predict, interact, N, shrunk.fit, AICc) {
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

    models = foreach(i=1:n, .packages=c('glmnet'), .errorhandling='remove') %dopar% {
        #Fit one location's model here
        loc = coords.unique[i,]
        gw = gweights[[i]]
		
		if (is.null(oracle)) {
	        m = gwglmnet.fit.inner(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, s=s, verbose=verbose, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, predict=predict, tuning=tuning, simulation=simulation, alpha=alpha, interact=interact, N=N, shrunk.fit=shrunk.fit, AICc=AICc)
		} else {
            m = gwlars.fit.oracle(x=x, y=y, family='gaussian', bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, AICc=AICc)
        }
        
        if (verbose) {
	        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", round(bw,3), "; s=", m[['s']], "; sigma2=", round(m[['sigma2']],3), "; nonzero=", paste(m[['nonzero']], collapse=","), "; weightsum=", round(m[['weightsum']],3), "; alpha=", round(m[['alpha']], 3), ".\n", sep=''))
        }
        return(m)
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
