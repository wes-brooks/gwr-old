gwglmnet.fit.fixedbwparallel = function(x, y, family, coords, bw, D=NULL, s=NULL, verbose=FALSE, mode.select, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE, alpha, simulation, tuning, predict, interact, N) {
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

        m = gwglmnet.fit.inner(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, s=s, verbose=verbose, mode.select=mode.select, gwr.weights=gw, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, predict=predict, tuning=tuning, simulation=simulation, alpha=alpha, interact=interact, N=N)
        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", bw, "; loss=", m[['loss.local']], ".\n", sep=''))

        return(m)
    }

    gwglmnet.object[['models']] = models
    gwglmnet.object[['coords']] = coords
    gwglmnet.object[['s.range']] = s

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
