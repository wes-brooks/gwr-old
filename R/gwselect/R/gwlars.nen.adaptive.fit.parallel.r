gwlars.nen.fit.parallel = function(x, y, coords, D, s, verbose, prior.weights, gweight, target, beta1, beta2, tol=1e-25, longlat=FALSE, adapt) {
    #Fit the GWLARS model (adaptive algorithm)
    coords.unique = unique(coords)

    gwlars.object = foreach(i=1:n, .packages='lars', .errorhandling='remove') %dopar% {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
        dist = D[i,]

        opt = optimize(gwlars.ssr, lower=beta1, upper=beta2, 
            maximum=FALSE, tol=target/1000, x=x, y=y, colocated=colocated, s=s,
            gweight=gweight, verbose=verbose, dist=dist, adapt=adapt,
            prior.weights=prior.weights, target=target, type=type)
        bandwidth = opt$minimum

        cat(paste("For i=", i, ", bw=", bandwidth, ", tolerance=", target/1000, ", miss=", sqrt(opt$objective), ".\n", sep=''))

        loow = gweight(D[i,-colocated], bandwidth)        
        prior.loow = prior.weights[-colocated]
        w <- prior.loow * loow
        reps = length(colocated)     

        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
           
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
        
        m <- ncol(xx)
        n <- nrow(xx)

        if (adapt==TRUE) {
            one <- rep(1, n)
            meanx <- drop(one %*% xx)/n
            x.centered <- scale(xx, meanx, FALSE)         # first subtracts mean
            normx <- sqrt(drop(one %*% (x.centered^2)))
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
            
            lm.step = try(lm(yy~xs, weights=w))  # mle fit on standardized
        
            if(class(lm.step) == "try-error") { 
                cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
                return(Inf)
            }
        
            beta.lm = lm.step$coeff[2:(m+1)]                    # mle except for intercept
            adapt.weight = abs(beta.lm)                        # weights for adaptive lasso
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(adapt.weight[k])) {
                    xs[,k] = xs[,k] * adapt.weight[k]
                } else {
                    xs[,k] = rep(0, dim(xs)[1])
                    adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
                }
            }
        } else {
            xs=xx
            meanx = rep(0, ncol(x))
            adapt.weight[k] = rep(1, ncol(x))
            normx[k] = rep(1, ncol(x))
        }

        xs.colocated = (x[colocated,] - meanx) * adapt.weight / normx
        model = lars(x=xs, y=yy, weights=w, type='lar')
        ll = model$lambda
        predictions = predict(model, newx=matrix(xs.colocated, nrow=reps, ncol=dim(xx)[2]), s=ll, type='fit')[['fit']]
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
        s.optimal = ll[which.min(cv.error)]
        
        if (verbose) { cat(paste(i, "\n", sep='')) }         
        list(model=model, cv.error=cv.error, s=s.optimal, index=i)
    }

    print("returning from gwlars.nen.adaptive.fit.parallel")
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
