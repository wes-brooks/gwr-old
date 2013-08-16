gwglmnet.sel = function(formula, data=list(), family, range=NULL, weights=NULL, coords, oracle=NULL, indx=NULL, adapt=FALSE, gweight=gwr.Gauss, s=NULL, method="dist", mode.select='AIC', verbose=FALSE, longlat=FALSE, tol=.Machine$double.eps^0.25, parallel=FALSE, alpha=1, precondition=FALSE, interact=FALSE, shrunk.fit=TRUE, AICc=FALSE) {
    if (!is.logical(adapt)) 
        stop("adapt must be logical")
    if (is.null(longlat) || !is.logical(longlat)) 
        longlat <- FALSE
    if (missing(coords)) 
        stop("Observation coordinates have to be given")
        
    if (!is.null(range)) {
        beta1 = min(range)
        beta2 = max(range)
    } else {
        if (method == "dist") {
            bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
            difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
            if (any(!is.finite(difmin))) 
                difmin[which(!is.finite(difmin))] <- 0
            beta1 <- difmin/1000
            beta2 <- difmin
        } else if (method == 'knn') {
            beta1 <- 0
            beta2 <- 1
        } else if (method == 'nen') {
            if (family=='binomial') {beta2 = sum(weights/(mean(y)*(1-mean(y))) * (y-mean(y))**2)}
            if (family=='poisson') {beta2 = sum(weights/(mean(y)) * (y-mean(y))**2)}
            beta1 = beta2/1000
        }
    }

    opt <- optimize(gwglmnet.cv.f, lower=beta1, upper=beta2, 
        maximum=FALSE, formula=formula, indx=indx, coords=coords, oracle=oracle, s=s, family=family, mode.select=mode.select,
        gweight=gweight, verbose=verbose, longlat=longlat, data=data, method=method, alpha=alpha, shrunk.fit=shrunk.fit,
        weights=weights, tol=tol, adapt=adapt, parallel=parallel, precondition=precondition, N=1, interact=interact, AICc=AICc)

    bdwt <- opt$minimum
    res <- bdwt
    res
}
