gwglmnet.nen.sel = function(formula, data=list(), coords, adapt=FALSE, gweight=gwr.Gauss, s, method="cv", verbose=FALSE, longlat=FALSE, family, weights, tol=.Machine$double.eps^0.25, type) {
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

    #Get the matrix of distances
    n = dim(coords)[1]
    if (longlat) {
        D = as.matrix(earth.dist(coords),n,n)
    } else {
        Xmat = matrix(rep(coords[,1], times=n), n, n)
        Ymat = matrix(rep(coords[,2], times=n), n, n)
        D = sqrt((Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2)
    }

    model = glm(formula=formula, data=data, family=family, weights=weights)
    difmin = sum(residuals(model, type=type)**2)

    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
        difmin[which(!is.finite(difmin))] <- 0
    beta1 <- difmin/1000
    beta2 <- difmin     

    opt <- optimize(gwglmnet.nen.cv.f, lower=beta1, upper=beta2, 
        maximum=FALSE, tol=tol, formula=formula, coords=coords, s=s,
        gweight=gweight, verbose=verbose, longlat=longlat, data=data, D=D,
        weights=weights, adapt=adapt, family=family, type=type)
    bdwt <- opt$minimum
    res <- bdwt

    res
}




gwglmnet.nen.cv.f = function(formula, data, bw, coords, gweight, verbose, adapt, longlat, s, family, weights, D=NULL, ...) {    
    #Generate the model with the given bandwidth:
    cat(paste('Beginning with bandwidth: ', bw, '\n', sep=''))
    gwglmnet.model = gwglmnet.nen(formula=formula, data=data, coords=coords, gweight=gweight, bw=bw, verbose=verbose, longlat=longlat, adapt=adapt, s=s, family=family, weights=weights, D=D )
    cv.error = sum(sapply(gwglmnet.model[['cv.error']], min))

    cat(paste('Bandwidth: ', bw, '. CV error: ', cv.error, '\n', sep=''))
    return(cv.error)
}


gwglmnet.nen <- function(formula, data, coords, gweight, bw, D=NULL, verbose=FALSE, longlat=FALSE, adapt=FALSE, s, family, weights=NULL) {
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
    n = dim(D)[1]

    if (!adapt) {
        res = gwglmnet.nen.fit(x, y, coords, bw, s, verbose, family, weights)
    }
    else {
        res = gwglmnet.nen.adaptive.fit(x, y, coords, bw, s, verbose, family, weights)
    }
    res[['data']] = data
    res[['response']] = as.character(formula[[2]])
    res
}


gwglmnet.nen.fit = function(x, y, coords, weight.matrix, s, verbose, family, prior.weights, gweight, bw) {
    #Fit the gwglmnet model (non-adaptive algorithm)
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    gwglmnet.object = list()
    cv.error = list()

    

    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])

        bandwidth = optimize(gwglmnet.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=tol, formula=formula, coords=coords, s=s,
            gweight=gweight, verbose=verbose, longlat=longlat, data=data, D=D,
            weights=weights, family=family, type=type)

        loow = weight.matrix[i,-colocated]
        prior.loow = prior.weights[-colocated]
        reps = length(colocated)        
        w <- prior.loow * loow
        
        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
        reps = length(colocated)        
        
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
        
        if (family=='binomial' && (abs(sum(yy*w)-sum(w))<1e-4 || sum(yy*w)<1e-4)) {
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yy*w), ", sum of weights=", sum(w), "\n", sep=''))
            model[[i]] = NULL
            cv.error[[i]] = 0
            s.optimal = c(s.optimal, max(s))
        } else {
            model[[i]] = glmnet(x=xx, y=cbind(1-yy, yy), weights=w, family=family, lambda=s)
            predictions = predict(model[[i]], newx=matrix(x[colocated,], nrow=reps, ncol=dim(xx)[2]), s=s, type='response')
            cv.error[[i]] = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(s)), nrow=reps, ncol=length(s))))
            s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
        }
        
        if (verbose) { cat(paste(i, "\n", sep='')) }
    }
    gwglmnet.object[['coef.scale']] = NULL
    gwglmnet.object[['model']] = model
    gwglmnet.object[['s']] = s.optimal
    gwglmnet.object[['mode']] = mode
    gwglmnet.object[['coords']] = coords.unique
    gwglmnet.object[['cv.error']] = cv.error
    gwglmnet.object[['s.range']] = s
    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}


gwglmnet.ssr <- function(bw, x, y, colocated, weight.matrix, s, verbose, family, prior.weights, gweight, D, type) {
    #Get the weights
    reps = length(colocated)
    loow = gweight(D, bw)[-colocated]
    w <- prior.weights[-colocated] * loow
    
    #Build the model
    xx = as.matrix(x[-colocated,])
    yy = as.matrix(y[-colocated])
    if (family=='binomial') {model = glmnet(x=xx, y=cbind(1-yy,yy), weights=w, family=family, lambda=s)}
    else {model = glmnet(x=xx, y=yy, weights=w, family=family, lambda=s)}

    #Find lambda to minimize CV error
    predictions = predict(model, newx=matrix(x[colocated,], nrow=reps, ncol=dim(xx)[2]), s=s, type='response')
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(s)), nrow=reps, ncol=length(s))))
    s.optimal = s[which.min(cv.error)]
    
    #Get the residuals at this choice of s (Poisson-specific for now):
    fitted = predict(model, newx=xx, s=s.optimal, type='response')
    pearson.resid = sum((yy - fitted)**2/fitted)

    pearson.resid
}


gwglmnet.nen.adaptive.fit = function(x, y, coords, weight.matrix, s, verbose, family, prior.weights) {
#Fit the gwglmnet model (adaptive algorithm)
    gwglmnet.object = list()
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    adapt.normx = list()
    adapt.scale = list()
    cv.error = list()
    coef.scale = list()
    glm.step = list()
    
    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
        loow = weight.matrix[i,-colocated]
        if (sum(loow)==0) { return(list(cv.error = Inf)) }      

        prior.loow = prior.weights[-colocated] 
        reps = length(colocated)        
        w <- prior.loow * loow
        
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])

        if (family=='binomial' && (abs(sum(yy*w)-sum(w))<1e-4 || sum(yy*w)<1e-4)) {            
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yy*w), ", sum of weights=", sum(w), "\n", sep=''))
            model[[i]] = NULL
            cv.error[[i]] = 0
            s.optimal = c(s.optimal, max(s))
        } else {
            m <- ncol(xx)
            n <- nrow(xx)
            one <- rep(1, n)
            meanx <- drop(one %*% xx)/n
            x.centered <- scale(xx, meanx, FALSE)         # first subtracts mean
            normx <- sqrt(drop(one %*% (x.centered^2)))
            adapt.normx[[i]] = normx
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
            
            out.glm = try(glm(yy~xs, family=family, weights=w))  # mle fit on standardized

            if(class(out.glm) == "try-error") { 
                cat(paste("Had to use the last glm for location ", i, "\n", sep=""))
                glm.step[[i]] = out.glm = glm.step[[i-1]]
            }
            else { glm.step[[i]] = out.glm }

            beta.glm = out.glm$coeff[2:(m+1)]                    # mle except for intercept
            adapt.weight = abs(beta.glm)                        # weights for adaptive lasso
            adapt.scale[[i]] = adapt.weight
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(adapt.weight[k])) {
                    xs[,k] = xs[,k] * adapt.weight[k]
                } else {
                    xs[,k] = rep(0, dim(xs)[1])
                    adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
                }
            }
            
            #Use the lars algorithm to fit the model
            coef.scale[[i]] = adapt.weight/normx
            names(coef.scale[[i]]) = sapply(strsplit(names(coef.scale[[i]]), 'xs'), function(x) {x[2]})
            
            if (sum(coef.scale[[i]]) <1e-10) {
                if (verbose) {cat(paste("opted for the intercept-only model at location: ", i, "\n", sep=""))}
                model[[i]] = NULL
                predictions = rep(coef(out.glm)[["(Intercept)"]], length(colocated))
                cv.error[[i]] = abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(s))))
                s.optimal = c(s.optimal, max(s))
            } else {
                if (family=='binomial') {model[[i]] = glmnet(x=xs, y=cbind(1-yy, yy), lambda=s, family=family, weights=w)}
                else {model[[i]] = glmnet(x=xs, y=yy, lambda=s, family=family, weights=w)}
                predictions = predict(model[[i]], newx=scale(matrix(x[colocated,], nrow=reps, ncol=dim(xx)[2]), center=meanx, scale=normx/adapt.weight), type='response', s=s)
                cv.error[[i]] = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(s)), nrow=reps, ncol=length(s))))
                s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
            }
        }

        if (verbose) { cat(paste(i, "\n", sep='')) }
    }
    gwglmnet.object[['coef.scale']] = coef.scale
    gwglmnet.object[['model']] = model
    gwglmnet.object[['s']] = s.optimal
    gwglmnet.object[['mode']] = mode
    gwglmnet.object[['coords']] = coords.unique
    gwglmnet.object[['cv.error']] = cv.error
    gwglmnet.object[['s.range']] = s
    class(gwglmnet.object) = 'gwglmnet.object'
    return(gwglmnet.object)
}
    