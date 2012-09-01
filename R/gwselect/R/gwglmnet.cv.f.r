gwglmnet.cv.f = function(formula, data, weights, family, bw, coords, gweight, verbose, adapt, longlat, s, tol, method, parallel, precondition) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family, weights=weights, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel, precondition=precondition)

    cv.error = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['cv.error']])}))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
