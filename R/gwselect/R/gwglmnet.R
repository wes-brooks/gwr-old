gwglmnet <-
function(formula, data, coords, gweight, bw, D=NULL, verbose=FALSE, longlat=FALSE, adapt=FALSE, s, family, weights=NULL, nearest.neighbors=FALSE) {
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

    if (is.null(D)) {
        #Get the matrix of distances
        n = dim(coords)[1]
        if (longlat) {
            D = as.matrix(earth.dist(coords),n,n)
        } else {
            Xmat = matrix(rep(coords[,1], times=n), n, n)
            Ymat = matrix(rep(coords[,2], times=n), n, n)
            D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
        }
    }    

    #Get the weight matrix
    if (!nearest.neighbors) {
        weight.matrix = gweight(D, bw)
    } else {
        n = dim(D)[1]
        bandwidths = sapply(1:n, function(x) {neighbor.weight(q=bw, D=D[x,], weight.function=gweight, verbose=verbose, tol=0.001)})
        weight.matrix = as.matrix(rbind(sapply(1:n, function(k) {gweight(as.vector(D[k,]), as.numeric(bandwidths[1,k]))})),n,n)
    }

    if (!adapt) {
        res = gwglmnet.fit(x, y, coords, weight.matrix, s, verbose, family, weights)
    }
    else {
        res = gwglmnet.adaptive.fit(x, y, coords, weight.matrix, s, verbose, family, weights)
    }
    res[['data']] = data
    res[['response']] = as.character(formula[[2]])
    res
}
