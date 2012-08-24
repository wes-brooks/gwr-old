gwlars.cv.f = function(formula, data, bw, coords, gweight, verbose, adapt, longlat, mode, s, tol, ...) {    
    #Generate the model with the given bandwidth:
    gwlars.model = gwlars(formula=formula, data=data, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s)
    cv.error = sum(sapply(gwlars.model[['cv.error']], min))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
