gwglmnet <- function(formula, data, family, weights=NULL, coords, fit.loc=NULL, indx=NULL, tuning=FALSE, predict=FALSE, simulation=FALSE, gweight, bw=NULL, mode.select='AIC', verbose=FALSE, longlat, tol, N=1, method, adapt=FALSE, s=NULL, parallel=FALSE, precondition=FALSE, D=NULL, interact=FALSE, alpha=1) {
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

    res = list()

    if (method=='dist') {
        weight.matrix = gweight(D, bw)
        if (parallel) {
            res[['model']] = gwglmnet.fit.fixedbwparallel(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, indx=indx, N=N, coords=coords, mode.select=mode.select, bw=bw, fit.loc=fit.loc, gwr.weights=weight.matrix, s=s, verbose=verbose, adapt=adapt, precondition=precondition, interact=interact, alpha=alpha)
        } else {
            res[['model']] = gwglmnet.fit.fixedbw(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, indx=indx, N=N, coords=coords, mode.select=mode.select, bw=bw, fit.loc=fit.loc, gwr.weights=weight.matrix, s=s, verbose=verbose, adapt=adapt, precondition=precondition, interact=interact, alpha=alpha)
        }
    } else {        
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        beta1 = difmin/300
        beta2 = difmin 

        if (method=='nen') {
            if (parallel) {
                res[['model']] = gwglmnet.fit.nenparallel(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, indx=indx, N=N, coords=coords, fit.loc=fit.loc, D=D, longlat=longlat, mode.select=mode.select, s=s, verbose=verbose, adapt=adapt, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol=tol, precondition=precondition, interact=interact, alpha=alpha)
            } else {
                res[['model']] = gwglmnet.fit.nen(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, indx=indx, coords=coords, fit.loc=fit.loc, N=N, D=D, longlat=longlat, mode.select=mode.select, s=s, verbose=verbose, adapt=adapt, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol=tol, precondition=precondition, interact=interact, alpha=alpha)
            }
        } else if (method=='knn') {
            if (parallel) {
                res[['model']] = gwglmnet.fit.knnparallel(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, indx=indx, coords=coords, fit.loc=fit.loc, N=N, D=D, longlat=longlat, mode.select=mode.select, s=s, verbose=verbose, adapt=adapt, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol=tol, precondition=precondition, interact=interact, alpha=alpha)
            } else {
                res[['model']] = gwglmnet.fit.knn(x=x, y=y, family=family, prior.weights=weights, tuning=tuning, predict=predict, simulation=simulation, indx=indx, coords=coords, fit.loc=fit.loc, N=N, D=D, longlat=longlat, mode.select=mode.select, s=s, verbose=verbose, adapt=adapt, target=bw, gweight=gweight, beta1=beta1, beta2=beta2, tol=tol, precondition=precondition, interact=interact, alpha=alpha)
            }
        }
    }

    if (!tuning) {
        res[['data']] = data
        res[['response']] = as.character(formula[[2]])
        res[['family']] = family
        res[['weights']] = weights
        res[['coords']] = coords
        res[['fit.locs']] = fit.loc
        res[['indx']] = indx
        res[['longlat']] = longlat
        res[['gweight']] = gweight
        res[['bw']] = bw
        res[['method']] = method
        res[['adapt']] = adapt
        res[['precondition']] = precondition
        res[['s']] = s
        res[['mode.select']] = mode.select
        res[['interact']] = interact
    }
    class(res) = "gwselect"

    res
}
