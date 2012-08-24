gwglmnet.cv.f <-
function(formula, data, bw, coords, gweight, verbose, adapt, longlat, s, family, weights, nn, D=NULL, ...) {    
    #Generate the model with the given bandwidth:
    cat(paste('Beginning with bandwidth: ', bw, '\n', sep=''))
    gwglmnet.model = gwglmnet(formula=formula, data=data, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s, family=family, weights=weights, nearest.neighbors=nn, D=D )
    cv.error = sum(sapply(gwglmnet.model[['cv.error']], min))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}
