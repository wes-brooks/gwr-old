neighbor.weight = function(q, D=NULL, this.coords=NULL, obs.coords=NULL, longlat=FALSE, weight.function, verbose=FALSE, tol=.Machine$double.eps^0.25) {
    if (is.null(D)) {
        bbox <- cbind(range(obs.coords[, 1]), range(obs.coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0

        this.coords = as.numeric(this.coords)
    
        if (!longlat) {
            nrow = dim(obs.coords)[1]
            D2 = (matrix(as.numeric(this.coords), nrow, 2, byrow=TRUE) - obs.coords)**2
            D2 = apply(D2, 1, sum)
            D = sqrt(D2)
        } else {
            n = dim(obs.coords)[1]
            D = sapply(1:dim(obs.coords)[1], function(x) {deg.dist(this.coords[1], this.coords[2], obs.coords[x,1], obs.coords[x,2])})
        }
    }

    beta1 <- min(D)
    beta2 <- 10*max(D)

    optimize(neighbor.diff, lower=beta1, upper=beta2, maximum=FALSE, tol=tol, 
                D=D, weight.function=weight.function, q.target=q, verbose=verbose)
}