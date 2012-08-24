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
    gwlars.object[['s.range']] = s
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}