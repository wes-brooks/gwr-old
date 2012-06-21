library(fossil)

gwlars.sel = function(formula, data=list(), coords, adapt=FALSE, gweight=gwr.Gauss, method="cv", verbose=TRUE, longlat=NULL, RMSE=FALSE, weights, tol=.Machine$double.eps^0.25) {
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
    m <- match(c("formula", "data", "weights"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    weights <- as.vector(model.extract(mf, "weights"))
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
    if (!adapt) {
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        beta1 <- difmin/1000
        beta2 <- difmin
        if (method == "cv") {
            opt <- optimize(gwlars.cv.f, lower = beta1, upper = beta2, 
                maximum = FALSE, y = y, x = x, coords = coords, 
                gweight = gweight, verbose = verbose, longlat = longlat, 
                RMSE = RMSE, weights = weights, tol = tol)
        }
        else {
            opt <- optimize(gwr.aic.f, lower = beta1, upper = beta2, 
                maximum = FALSE, y = y, x = x, coords = coords, 
                gweight = gweight, verbose = verbose, longlat = longlat, 
                tol = tol)
        }
        bdwt <- opt$minimum
        res <- bdwt
    }
    else {
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        beta1 <- difmin/1000
        beta2 <- difmin
        if (method == "cv") {
            opt <- optimize(gw.adaptive.lars.cv.f, lower = beta1, upper = beta2, 
                maximum = FALSE, y = y, x = x, coords = coords, 
                gweight = gweight, verbose = verbose, longlat = longlat, 
                RMSE = RMSE, weights = weights, tol = tol)
        }
        else {
            opt <- optimize(gwr.aic.f, lower = beta1, upper = beta2, 
                maximum = FALSE, y = y, x = x, coords = coords, 
                gweight = gweight, verbose = verbose, longlat = longlat, 
                tol = tol)
        }
        bdwt <- opt$minimum
        res <- bdwt
    }
    res
}




gwlars.cv.f = function(bw, y, x, coords, gweight, verbose, longlat, tol, ...) {    
    cat(paste('Entering gwlars.cv.f with bandwidth ', bw, '\n', sep=''))
    #Compute the matrix of distances (in kilometers)
    n = dim(coords)[1]
    D = as.matrix(earth.dist(coords),n,n)
    w = bisquare(D, bw=bw)
    
    #Do the same for just the unique locations
    coords.unique = unique(coords)

    gwlars = list()
    lambda = seq(0, 5, length.out=2000)
    l = vector()
    error.cv = 0

    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
        loow = w[i,-colocated]
    
        w.sqrt <- diag(sqrt(loow))
        gwlars[[i]] = lars(x=w.sqrt %*% as.matrix(x[-colocated,]), y=w.sqrt %*% as.matrix(y[-colocated]))
        
        errs = colSums(abs(predict(gwlars[[i]], newx=x[colocated,], s=lambda, type='fit', mode='lambda')[['fit']] - y[colocated]))
        error.cv = error.cv + min(errs)
        l = c(l, which.min(errs)/1000)
        cat(paste(i, '\n', sep=''))
    }

    cat(paste('Bandwidth: ', bw, '. CV error: ', error.cv, '\n', sep=''))
    return(error.cv)
}


gw.adaptive.lars.cv.f = function(bw, y, x, coords, gweight, verbose, longlat, tol, ...) {    
    cat(paste('Entering gwadalars.cv.f with bandwidth ', bw, '\n', sep=''))
    #Compute the matrix of distances (in kilometers)
    n = dim(coords)[1]
    D = as.matrix(earth.dist(coords),n,n)
    w = bisquare(D, bw=bw)
    
    #Do the same for just the unique locations
    coords.unique = unique(coords)

    gwlars = list()
    lambda = seq(0, 5, length.out=2000)
    l = vector()
    error.cv = 0

    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
        loow = w[i,-colocated]
    
        #Apply the geographic weights:
        w.sqrt <- diag(sqrt(loow))
        x.weighted = w.sqrt %*% as.matrix(x[-colocated,])
        y.weighted = w.sqrt %*% as.matrix(y[-colocated])

        #Generate the adaptive lasso weights
        m<-ncol(x.weighted)
        n<-nrow(x.weighted)
        one <- rep(1, n)
        meanx <- drop(one %*% x.weighted)/n
        x.centered <- scale(x.weighted, meanx, FALSE)         # first subtracts mean
        normx <- sqrt(drop(one %*% (x.centered^2)))
        names(normx) <- NULL
        
        xs = x.centered
        for (k in 1:dim(x.centered)[2]) {
            if (normx[k]!=0) {
                xs[,k] = xs[,k] / normx[k]
            } else {
                xs[,k] = rep(0, dim(xs)[1])
            }
        }
    
        out.ls = lm(y.weighted~xs)                      # ols fit on standardized
        beta.ols = out.ls$coeff[2:(m+1)]       # ols except for intercept
        ada.weight = abs(beta.ols)                      # weights for adaptive lasso
        for (k in 1:dim(x.centered)[2]) {
            if (!is.na(ada.weight[k])) {
                xs[,k] = xs[,k] * ada.weight[k]
            } else {
                xs[,k] = rep(0, dim(xs)[1])
            }
        }
        
        #Use the lars algorithm to fit the model
        gwlars[[i]] = lars(x=xs, y=y.weighted, normalize=FALSE)

        errs = colSums(abs(predict(gwlars[[i]], newx=scale(x[colocated,], center=meanx, scale=normx/ada.weight), s=lambda, type='fit', mode='lambda')[['fit']] - y[colocated,]))
        error.cv = error.cv + min(errs)
        l = c(l, which.min(errs)/1000)
        cat(paste(i, '\n', sep=''))
    }

    cat(paste('Bandwidth: ', bw, '. CV error: ', error.cv, '\n', sep=''))
    return(error.cv)
}