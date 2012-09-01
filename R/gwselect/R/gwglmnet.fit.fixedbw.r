gwglmnet.fit.fixedbw = function(x, y, family, coords, bw, D=NULL, s=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, precondition=FALSE) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwglmnet.object = list()
    models = list()

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(D, bw)    
    }       

    for (i in 1:n) {
        #Fit one location's model here
        loc = coords.unique[i,]
        models[[i]] = gwglmnet.fit.inner(x=x, y=y, family=family, bw=bw, coords=coords, loc=loc, dist=D[i,], s=s, verbose=verbose, gwr.weights=gwr.weights[i,], prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition)
    }

    gwglmnet.object[['models']] = models
    gwglmnet.object[['coords']] = coords
    gwglmnet.object[['s.range']] = s

    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
