gwglmnet.cv.f = function(formula, data, weights, indx, family, bw, coords, gweight, mode, mode.select, verbose, adapt, longlat, s, tol, method, N, parallel, precondition, shrink) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family, weights=weights, indx=indx, coords=coords, gweight=gweight, bw=bw, mode=mode, N=N, mode.select=mode.select, verbose=verbose, longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel, precondition=precondition, shrink=shrink)

    cv.error = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['loss.local']])}))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
