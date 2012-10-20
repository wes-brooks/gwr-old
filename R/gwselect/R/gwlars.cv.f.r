gwlars.cv.f = function(formula, data, weights, bw, N=N, coords, gweight, verbose, adapt, longlat, mode, s, mode.select, shrink, tol, method, parallel, precondition) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwlars.model = gwlars(formula=formula, data=data, weights=weights, coords=coords, longlat=longlat, gweight=gweight, N=N, bw=bw, adapt=adapt, mode=mode, s=s, mode.select=mode.select, shrink=shrink, method=method, parallel=parallel, precondition=precondition, verbose=verbose)

    loss = sum(sapply(gwlars.model[['model']][['models']], function(x) {x[['loss.local']]}))

    cat(paste('Bandwidth: ', bw, '. Loss: ', loss, '\n', sep=''))
    return(loss)
}
