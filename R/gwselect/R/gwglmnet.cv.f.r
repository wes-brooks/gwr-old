gwglmnet.cv.f = function(formula, data, weights, indx, family, bw, coords, gweight, oracle, mode.select, verbose, adapt, longlat, s, tol, method, N, parallel, precondition, interact, alpha, shrunk.fit, AICc) {    
    #Generate the model with the given bandwidth:
    cat(paste("starting bw:", round(bw, 3), '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, family=family, weights=weights, tuning=TRUE, indx=indx, coords=coords, gweight=gweight, oracle=oracle, bw=bw, N=N, mode.select=mode.select, verbose=verbose, longlat=longlat, adapt=adapt, s=s, method=method, parallel=parallel, precondition=precondition, interact=interact, alpha=alpha, shrunk.fit=shrunk.fit, AICc=AICc)

    if (AICc) {
        trH = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['loss.local']]})) 

        #Local sigma^2:
        #loss = nrow(data) * (mean(sapply(gwglmnet.model[['model']][['models']], function(x) {log(x[['sigma2']])})) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))

        #Global sigma^2:
        #if (family=='gaussian') {
        #	#Use the AICc
        	loss = nrow(data) * (log(mean(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['ssr.local']]}))) + 1 + (2*(trH+1))/(nrow(data)-trH-2) + log(2*pi))
        #}
        #else if (family %in% c('binomial', 'poisson')) {
        #	#Use GCV (OSullivan et al. 1986)
        #	loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['ssr.local']]})) / (nrow(data)-trH)**2
        #}
        
        #"Simplistic" BIC - based on eq4.22 from the Fotheringham et al. book:
        #loss = nrow(data) * (log(mean(sapply(gwglmnet.model[['model']][['models']], function(x) {x[['ssr.local']]}))) + 1 + log(2*pi)) + trH * log(nrow(data))/2
    }
    else { loss = sum(sapply(gwglmnet.model[['model']][['models']], function(x) {min(x[['loss.local']])})) }

    cat(paste('Bandwidth: ', round(bw, 3), '. Loss: ', signif(loss, 5), '\n', sep=''))
    return(loss)
}
