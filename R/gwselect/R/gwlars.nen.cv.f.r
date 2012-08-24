gwlars.nen.cv.f = function(formula, data, bw, coords, gweight, verbose, adapt, longlat, mode, s, tol, parallel=FALSE, ...) {    
    #Generate the model with the given bandwidth:
    gwlars.model = gwlars.nen(formula=formula, data=data, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s, parallel=parallel)
    cv.error = sum(sapply(gwlars.model[['cv.error']], min))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
