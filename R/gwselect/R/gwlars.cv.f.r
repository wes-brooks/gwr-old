gwlars.cv.f = function(formula, data, weights, bw, N=N, coords, fit.loc, indx, gweight, verbose, adapt, longlat, mode, s, mode.select, tol, method, parallel, precondition, oracle) {    
    #Generate the model with the given bandwidth:
    cat(paste("preparing for bw:", bw, '\n', sep=''))
    gwlars.model = gwlars(formula=formula, data=data, weights=weights, coords=coords, fit.loc=fit.loc, longlat=longlat, gweight=gweight, N=N, bw=bw, indx=indx, tuning=TRUE, predict=FALSE, adapt=adapt, mode=mode, s=s, mode.select=mode.select, method=method, parallel=parallel, precondition=precondition, oracle=oracle, verbose=verbose)

    loss = sum(sapply(gwlars.model[['model']][['models']], function(x) {x[['loss.local']]}))

    cat(paste('Bandwidth: ', bw, '. Loss: ', loss, '\n', sep=''))
    return(loss)
}
