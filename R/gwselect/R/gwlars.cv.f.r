gwlars.cv.f = function(formula, data, weights, bw, coords, gweight, verbose, adapt, longlat, mode, s, tol, method, parallel) {    
    #Generate the model with the given bandwidth:
    gwlars.model = gwlars(formula=formula, data=data, weights=weights, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel)
    
    print(gwlars.model[['models']])

    cv.error = sum(sapply(gwlars.model[['models']], function(x) {min(x[['cv.error']])}))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
