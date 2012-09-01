gwlars.cv.f = function(formula, data, weights, bw, coords, gweight, verbose, adapt, longlat, mode, s, tol, method, parallel, precondition) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwlars.model = gwlars(formula=formula, data=data, weights=weights, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, mode=mode, s=s, method=method, parallel=parallel, precondition=precondition)

    cv.error = sum(sapply(gwlars.model[['model']][['models']], function(x) {min(x[['cv.error']])}))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
