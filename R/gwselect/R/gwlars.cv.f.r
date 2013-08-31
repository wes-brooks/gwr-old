gwlars.cv.f = function(formula, data, weights, bw, N=N, coords, fit.loc, indx, gweight, verbose, adapt, interact, longlat, s, mode.select, tol, method, parallel, precondition, oracle, shrunk.fit, AICc) {    
    #Generate the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    gwlars.model = gwlars(formula=formula, data=data, weights=weights, coords=coords, fit.loc=fit.loc, longlat=longlat, gweight=gweight, N=N, bw=bw, indx=indx, tuning=TRUE, predict=FALSE, adapt=adapt, s=s, mode.select=mode.select, method=method, parallel=parallel, precondition=precondition, oracle=oracle, verbose=verbose, interact=interact, shrunk.fit=shrunk.fit, AICc=AICc)

    if (AICc) {
        trH = sum(sapply(gwlars.model[['model']][['models']], function(x) {x[['loss.local']]})) 
        #print(log(mean(sapply(gwlars.model[['model']][['models']], function(x) {x[['sigma2']]}))))
        #print(nrow(data))
        #print(trH)
        #print((2*(trH+1))/(nrow(data)-trH-2))
        
        #Local sigma^2:
        #loss = nrow(data) * (mean(sapply(gwlars.model[['model']][['models']], function(x) {log(x[['sigma2']])})) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))

        #Global sigma^2:
        loss = nrow(data) * (log(mean(sapply(gwlars.model[['model']][['models']], function(x) {x[['ssr.local']]}))) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))
    }
    else {loss = sum(sapply(gwlars.model[['model']][['models']], function(x) {x[['loss.local']]}))}

    cat(paste('Bandwidth: ', round(bw, 3), '. Loss: ', round(loss, 3), '\n', sep=''))
    return(loss)
}
