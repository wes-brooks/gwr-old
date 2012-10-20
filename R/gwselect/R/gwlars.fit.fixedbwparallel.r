gwlars.fit.fixedbwparallel = function(x, y, coords, bw, D=NULL, s=NULL, mode.select, shrink, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, mode, precondition=FALSE) {
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwlars.object = list()
    models = list()

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(D, bw)    
    }       

    models = foreach(i=1:n, .packages=c('lars'), .errorhandling='remove') %dopar% {
        #Fit one location's model here
        loc = coords.unique[i,]
        m = gwlars.fit.inner(x=x, y=y, bw=bw, coords=coords, loc=loc, dist=D[i,], N=N, s=s, mode.select=mode.select, shrink=shrink, verbose=verbose, gwr.weights=gwr.weights[i,], prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode=mode, precondition=precondition)
        cat(paste("For i=", i, ", bw=", bw, ", loss=", paste(m[['loss']], collapse=","), "(min=", m[['loss']][m[['s']]], ").\n", sep=''))
        
        return(m)
    }

    gwlars.object[['models']] = models
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['s.range']] = s

    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
