gwlars.fit.fixedbw = function(x, y, coords, bw, D=NULL, s=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, mode) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwlars.object = list()
    models = list()

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(D, bw)    
    }       

    for (i in 1:n) {
        #Fit one location's model here
        loc = coords.unique[i,]
        models[[i]] = gwlars.fit.inner(x=x, y=y, coords=coords, loc=loc, dist=NULL, s=s, verbose=verbose, gwr.weights=gwr.weights[i,], prior.weights=prior.weights, gweight=NULL, adapt=adapt, mode=mode)
    }

    gwlars.object[['models']] = models
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['s.range']] = s

    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
