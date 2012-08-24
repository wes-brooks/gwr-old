gwglmnet.nen.cv.f <-
function(formula, data, bw, coords, gweight, verbose, adapt, longlat, s=NULL, beta1, beta2, family, weights=NULL, D=NULL, tolerance=.Machine$double.eps^0.25, type='pearson', parallel=FALSE, ...) {    
    #Generate the model with the given bandwidth:
    cat(paste('Beginning with target SSR: ', bw, ', tolerance: ', tolerance, '\n', sep=''))
    gwglmnet.model = gwglmnet.nen(formula=formula, data=data, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s, family=family, weights=weights, D=D, tol=tolerance, beta1=beta1, beta2=beta2, type, parallel=parallel)

    cv.error = sum(sapply(gwglmnet.model[['model']], function(x) min(x[['cv.error']])))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
