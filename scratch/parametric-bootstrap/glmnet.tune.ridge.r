gwglmnet.tune.ridge = function(models, x, y, family, S, coords, fit.loc=NULL, oracle, bw, D=NULL, verbose=FALSE, gwr.weights=NULL, indx=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, precondition=FALSE, interact, N) {
    if (!is.null(fit.loc)) {
        coords.unique = fit.loc
    } else {
        coords.unique = unique(coords)
    }
    n = dim(coords.unique)[1]

    gwglmnet.object = list()
    resampled = vector()

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(D, bw)    
    }       

    gweights = list()
    for (j in 1:nrow(gwr.weights)) {
        gweights[[j]] = as.vector(gwr.weights[j,])
    }      

    quant = vector()

    for (i in 1:n) {
        #Fit one location's model here
        loc = coords.unique[i,]
        gw = gweights[[i]]

		if (is.null(oracle)) {
	        resampled = rbind(resampled, gwglmnet.fit.ridge(x=x, y=y, family=family, S=S, bw=bw, coords=coords, loc=loc, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, precondition=precondition, interact=interact))
		} else {
            resampled = rbind(resampled, gwselect.fit.ridge(x=x, y=y, family=family, S=S, bw=bw, coords=coords, loc=loc, indx=indx, oracle=oracle[[i]], N=N, mode.select=mode.select, tuning=tuning, predict=predict, simulation=simulation, verbose=verbose, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, AICc=AICc))
        }
        
        ee = ecdf(resampled[nrow(resampled),])
        quant = c(quant, ee(y[i]))
        if (verbose) {            
	        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); quantile=", tail(quant, 1), "; resampled=", paste(round(resampled[nrow(resampled),],3), collapse=" "), "; truth=", round(y[i],3), ".\n", sep=''))
		}
    }
	
    return(resampled)
}
