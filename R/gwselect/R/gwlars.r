gwlars <- function(formula, data, weights=NULL, coords, gweight, bw=NULL, verbose=FALSE, longlat, tol, method, adapt=FALSE, s=NULL, mode='lambda', parallel=FALSE) {
    if (!is.logical(adapt)) 
        stop("adapt must be logical")
    if (is(data, "Spatial")) {
        if (!missing(coords)) 
            warning("data is Spatial* object, ignoring coords argument")
        coords <- coordinates(data)
        if ((is.null(longlat) || !is.logical(longlat)) && !is.na(is.projected(data)) && 
            !is.projected(data)) {
            longlat <- TRUE
        }
        else longlat <- FALSE
        data <- as(data, "data.frame")
    }
    if (is.null(longlat) || !is.logical(longlat)) 
        longlat <- FALSE
    if (missing(coords)) 
        stop("Observation coordinates have to be given")
    mf <- match.call(expand.dots = FALSE)
    #m <- match(c("formula", "data", "weights"), names(mf), 0)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    #weights <- as.vector(model.extract(mf, "weights"))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) 
        weights <- rep(as.numeric(1), dp.n)
    if (any(is.na(weights))) 
        stop("NAs in weights")
    if (any(weights < 0)) 
        stop("negative weights")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)

    #Get the matrices of distances and weights
    n = dim(coords)[1]
    if (longlat) {
        D = as.matrix(earth.dist(coords),n,n)
    } else {
        Xmat = matrix(rep(coords[,1], times=n), n, n)
        Ymat = matrix(rep(coords[,2], times=n), n, n)
        D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
    }

    res = list()

    if (method=='nen') {
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        beta1 = difmin/1000
        beta2 <- difmin/3   

        if (parallel) {
            print("going in parallel")
            res[['models']] = gwlars.fit.nenparallel(x=x, y=y, prior.weights=weights, coords=coords, D=D, longlat=longlat, s=s, mode=mode, verbose=verbose, adapt=adapt, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol=tol)
        } else {
            res[['models']] = gwlars.fit.nen(x=x, y=y, prior.weights=weights, coords=coords, D=D, longlat=longlat, s=s, mode=mode, verbose=verbose, adapt=adapt, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol=tol)
        }
    } else if (method=='distance') {
        weight.matrix = gweight(D, bw)
        res[['models']] = gwlars.fit.fixedbw(x=x, y=y, weights=weights, coords=coords, weight.matrix=weight.matrix, s=s, mode=mode, verbose=verbose, adapt=adapt)
    } else if (method=='knn') {
        res[['models']] = gwlars.fit.knn(x=x, y=y, weights=weights, coords=coords, weight.matrix=weight.matrix, s=s, mode=mode, verbose=verbose, adapt=adapt)
    }


    res[['data']] = data
    res[['response']] = as.character(formula[[2]])
    res
}
