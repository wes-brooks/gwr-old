gwglmnet.cv.f = function(formula, data, weights, indx, family, bw, coords, gweight, mode.select, verbose, adapt, longlat, s, tol, method, N, parallel, precondition, interact, alpha, shrunk.fit) {    
    #Generate the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family, weights=weights, tuning=TRUE, indx=indx, coords=coords, gweight=gweight, bw=bw, N=N, mode.select=mode.select, verbose=verbose, longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel, precondition=precondition, interact=interact, alpha=alpha, shrunk.fit=shrunk.fit)

    loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['loss.local']])}))

    cat(paste('Bandwidth: ', round(bw, 3), '. Loss: ', round(loss, 3), '\n', sep=''))
    return(loss)
}
