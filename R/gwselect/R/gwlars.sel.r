gwlars.sel = function(formula, data=list(), coords, indx=NULL, fit.loc=NULL, range=NULL, adapt=FALSE, gweight=gwr.Gauss, mode, s, N=1, mode.select="CV", shrink=TRUE, method="dist", verbose=FALSE, longlat=FALSE, weights=NULL, tol=.Machine$double.eps^0.25, parallel=FALSE, precondition=FALSE) {
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
            beta2 = sum(weights * (y-mean(y))**2)
            beta1 = beta2/1000
        }
    }

    opt <- optimize(gwlars.cv.f, lower=beta1, upper=beta2, 
        maximum=FALSE, formula=formula, coords=coords, indx=indx, s=s, N=N, mode=mode, mode.select=mode.select,
        gweight=gweight, verbose=verbose, longlat=longlat, data=data, method=method, shrink=shrink,
        weights=weights, tol=tol, adapt=adapt, parallel=parallel, precondition=precondition)

    bdwt <- opt$minimum
    res <- bdwt
    res
}
