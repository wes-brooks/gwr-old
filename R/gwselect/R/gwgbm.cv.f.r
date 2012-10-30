gwgbm.cv.f = function(formula, data, weights, family, bw, coords, gweight, verbose, adapt, longlat, s, tol, method, parallel, precondition) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwgbm.model = gwgbm(formula=formula, data=data, family=family, weights=weights, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, method=method, parallel=parallel, precondition=precondition)

    loss = sum(sapply(gwgbm.model[['model']][['models']], function(x) {min(x[['loss.local']])}))

    cat(paste('Bandwidth: ', bw, '. Loss: ', loss, '\n', sep=''))
    return(loss)
}
