gwglmnet.ridge <- function(formula, data, family, weights=NULL, coords, S=1, fit.loc=NULL, indx=NULL, oracle=NULL, gweight, bw=NULL, verbose=FALSE, longlat, tol, N=1, parallel=FALSE, precondition=FALSE, D=NULL, interact=FALSE) {
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
    if (is.null(D)) {
        n = dim(coords)[1]
        if (longlat) {
            D = as.matrix(earth.dist(coords),n,n)
        } else {
            Xmat = matrix(rep(coords[,1], times=n), n, n)
            Ymat = matrix(rep(coords[,2], times=n), n, n)
            D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
        }
    }

    weight.matrix = gweight(D, bw)

    res = gwglmnet.tune.ridge(x=x, y=y, family=family, S=S, prior.weights=weights, indx=indx, N=N, coords=coords, oracle=oracle, bw=bw, fit.loc=fit.loc, gwr.weights=weight.matrix, verbose=verbose, precondition=precondition, interact=interact)
    res
}
