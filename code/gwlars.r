gwlars <- function(formula, data, coords, gweight, bw, verbose=FALSE, longlat, tol, adapt=FALSE, s, mode='lambda') {
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

    #Get the matrices of distances and weights
    n = dim(coords)[1]
    if (longlat) {
        D = as.matrix(earth.dist(coords),n,n)
    } else {
        Xmat = matrix(rep(coords[,1], times=n), n, n)
        Ymat = matrix(rep(coords[,2], times=n), n, n)
        D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
    }
    weight.matrix = gweight(D, bw)

    if (!adapt) {
        res = gwlars.fit(x, y, coords, weight.matrix, s, mode, verbose)
    }
    else {
        res = gwlars.adaptive.fit(x, y, coords, weight.matrix, s, mode, verbose)
    }
    res
}


gwlars.fit = function(x, y, coords, weight.matrix, s, mode, verbose) {
#Fit the GWLARS model (non-adaptive algorithm)
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    gwlars.object = list()
    cv.error = list()
    coef.scale = list()
    
    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords$x==coords.unique$x[i] & coords$y==coords.unique$y[i])
        loow = weight.matrix[i,-colocated]
        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
        reps = length(colocated)        
        w.sqrt <- diag(rep(sqrt(loow), reps))
        
        coef.scale[[i]] = rep(1, dim(x)[2])
        model[[i]] = lars(x=w.sqrt %*% as.matrix(x[-colocated,]), y=w.sqrt %*% as.matrix(y[-colocated]))
        predictions = predict(model[[i]], newx=matrix(x[colocated,], nrow=reps, ncol=dim(x)[2]), s=s, type='fit', mode=mode)[['fit']]
        cv.error[[i]] = colSums(abs(matrix(predictions - as.matrix(y[colocated]), nrow=reps, ncol=length(s))))
        s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
        

        if (verbose) { cat(paste(i, "\n", sep='')) }
    }
    gwlars.object[['coef.scale']] = coef.scale
    gwlars.object[['model']] = model
    gwlars.object[['s']] = s.optimal
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['cv.error']] = cv.error
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}


gwlars.adaptive.fit = function(x, y, coords, weight.matrix, s, mode, verbose) {
#Fit the GWLARS model (adaptive algorithm)
    gwlars.object = list()
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    adalars.normx = list()
    adalars.scale = list()
    cv.error = list()
    coef.scale = list()
    
    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords$x==coords.unique$x[i] & coords$y==coords.unique$y[i])
        loow = weight.matrix[i,-colocated]
        if (sum(loow)==0) { return(list(cv.error = Inf)) }            
        reps = length(colocated)        
        w.sqrt <- diag(rep(sqrt(loow), reps))
        
        x.weighted = w.sqrt %*% x[-colocated,]
        y.weighted = w.sqrt %*% as.matrix(y[-colocated])
        m <- ncol(x.weighted)
        n <- nrow(x.weighted)
        one <- rep(1, n)
        meanx <- drop(one %*% x.weighted)/n
        x.centered <- scale(x.weighted, meanx, FALSE)         # first subtracts mean
        normx <- sqrt(drop(one %*% (x.centered^2)))
        adalars.normx[[i]] = normx
        names(normx) <- NULL
        xs = x.centered
        for (k in 1:dim(x.centered)[2]) {
            if (normx[k]!=0) {
                xs[,k] = xs[,k] / normx[k]
            } else {
                xs[,k] = rep(0, dim(xs)[1])
                normx[k] = Inf #This should allow the lambda-finding step to work.
            }
        }
        
        out.ls = lm(y.weighted~xs)                      # ols fit on standardized
        beta.ols = out.ls$coeff[2:(m+1)]       # ols except for intercept
        ada.weight = abs(beta.ols)                      # weights for adaptive lasso
        adalars.scale[[i]] = ada.weight
        for (k in 1:dim(x.centered)[2]) {
            if (!is.na(ada.weight[k])) {
                xs[,k] = xs[,k] * ada.weight[k]
            } else {
                xs[,k] = rep(0, dim(xs)[1])
                ada.weight[k] = 0 #This should allow the lambda-finding step to work.
            }
        }
        
        #Use the lars algorithm to fit the model
        coef.scale[[i]] = ada.weight/normx
        names(coef.scale[[i]]) = sapply(strsplit(names(coef.scale[[i]]), 'xs'), function(x) {x[2]})
        model[[i]] = lars(x=xs, y=y.weighted, normalize=FALSE)
        predictions = predict(model[[i]], newx=scale(matrix(x[colocated,], nrow=reps, ncol=dim(x)[2]), center=meanx, scale=normx/ada.weight), type='fit', mode=mode, s=s)[['fit']]
        cv.error[[i]] = colSums(abs(matrix(predictions - matrix(y[colocated]), nrow=reps, ncol=length(s))))
        s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])

        if (verbose) { cat(paste(i, "\n", sep='')) }
    }
    gwlars.object[['coef.scale']] = coef.scale
    gwlars.object[['model']] = model
    gwlars.object[['s']] = s.optimal
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['cv.error']] = cv.error
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
    