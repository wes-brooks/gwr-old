gwglmnet.fit.fixedbw = function(x, y, family, coords, bw, D=NULL, s=NULL, verbose=FALSE, mode.select, gwr.weights=NULL, indx=NULL, prior.weights=NULL, tuning=FALSE, predict=FALSE, fit.loc=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE, alpha, simulation, interact, N) {
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

        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, s=s, mode.select=mode.select, verbose=verbose, prior.weights=prior.weights, predict=predict, tuning=tuning, simulation=simulation, gwr.weights=gweights[[i]], adapt=adapt, precondition=precondition, alpha=alpha, interact=interact, N=N)
        cat(paste("For i=", i, "; location=(", paste(round(loc,3), collapse=","), "); bw=", bw, "; loss=", models[[i]][['loss.local']], ".\n", sep=''))
    }

    gwglmnet.object[['models']] = models
    gwglmnet.object[['coords']] = coords
    gwglmnet.object[['s.range']] = s

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
