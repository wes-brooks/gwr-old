gwglmnet.cv.f = function(formula, data, weights, indx, family, bw, coords, gweight, mode.select, verbose, adapt, longlat, s, tol, method, N, parallel, precondition, interact) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family, weights=weights, tuning=TRUE, indx=indx, coords=coords, gweight=gweight, bw=bw, N=N, mode.select=mode.select, verbose=verbose, longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel, precondition=precondition, interact=interact)

    loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['loss.local']])}))

    cat(paste('Bandwidth: ', bw, '. Loss: ', loss, '\n', sep=''))
    return(loss)
}
