gwlars.fit.parallel = function(x, y, weights, coords, weight.matrix, s, mode, verbose, adapt) {
    coords.unique = unique(coords)
    
    gwlars.object = foreach(i=1:n, .packages='lars', .errorhandling='remove') %dopar% {
        colocated = which(coords$x==coords.unique$x[i] & coords$y==coords.unique$y[i])
        loow = weight.matrix[i,-colocated]
        w <- prior.loow[-colocated] * loow
        if (sum(loow)==0) { return(list(cv.error = Inf)) }            
        reps = length(colocated)        
        sqrt.w <- diag(sqrt(w))

        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])

        m <- ncol(x.weighted)
        n <- nrow(x.weighted)

        if (adapt==TRUE) {
            one <- rep(1, n)
            meanx <- drop(one %*% xx)/n
            x.centered <- scale(xx, meanx, FALSE)         # first subtracts mean
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
            
            out.ls = lm(yy~xx, weights=w)                      # ols fit on standardized
            beta.ols = out.ls$coeff[2:(m+1)]                # ols except for intercept
            ada.weight = abs(beta.ols)                      # weights for adaptive lasso
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(ada.weight[k])) {
                    xs[,k] = xs[,k] * ada.weight[k]
                } else {
                    xs[,k] = rep(0, dim(xs)[1])
                    ada.weight[k] = 0 #This should allow the lambda-finding step to work.
                }
            }
            predx = scale(matrix(x[colocated,], nrow=reps, ncol=dim(x)[2]), center=meanx, scale=normx/ada.weight)            
        } else {
            xs = x.weighted
            normx = rep(1, ncol(xs))
            ada.weight = rep(1, ncol(xs))
            predx = as.matrix(x[colocated,])
        }
        
        #Use the lars algorithm to fit the model
        coef.scale = ada.weight/normx
        names(coef.scale) = sapply(strsplit(names(coef.scale), 'xs'), function(x) {x[2]})
        model = lars(x=xs, y=y.weighted, normalize=FALSE, type='lar')

        if (is.null(s)) {
            s = model$lambda
            mode='lambda'
        }
        predictions = predict(model, newx=predx, type='fit', mode=mode, s=s)[['fit']]
        cv.error = colSums(abs(matrix(predictions - matrix(y[colocated]), nrow=reps, ncol=length(s))))
        s.optimal = c(s.optimal, s[which.min(cv.error)])

        if (verbose) { cat(paste(i, "\n", sep='')) }         
        list(model=model, cv.error=cv.error, s=s.optimal, index=i)
    }

    gwlars.object[['coords']] = coords
    gwlars.object[['mode']] = mode
    gwlars.object[['s.range']] = s
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}
