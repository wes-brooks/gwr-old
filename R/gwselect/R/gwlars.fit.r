gwlars.fit = function(x, y, coords, weight.matrix, s, mode, verbose) {
#Fit the GWLARS model (non-adaptive algorithm)
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    gwlars.object = list()
    cv.error = list()
    
    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords$x==coords.unique$x[i] & coords$y==coords.unique$y[i])
        loow = weight.matrix[i,-colocated]
        if (sum(loow)==0) { return(list(cv.error = Inf)) }   
        reps = length(colocated)        
        w.sqrt <- diag(rep(sqrt(loow), reps))
        
        model[[i]] = lars(x=w.sqrt %*% as.matrix(x[-colocated,]), y=w.sqrt %*% as.matrix(y[-colocated]))
        predictions = predict(model[[i]], newx=matrix(x[colocated,], nrow=reps, ncol=dim(x)[2]), s=s, type='fit', mode=mode)[['fit']]
        cv.error[[i]] = colSums(abs(matrix(predictions - as.matrix(y[colocated]), nrow=reps, ncol=length(s))))
        s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
        

        if (verbose) { cat(paste(i, "\n", sep='')) }
    }
    gwlars.object[['coef.scale']] = NULL
    gwlars.object[['model']] = model
    gwlars.object[['s']] = s.optimal
    gwlars.object[['mode']] = mode
    gwlars.object[['coords']] = coords
    gwlars.object[['cv.error']] = cv.error
    gwlars.object[['s.range']] = s
    class(gwlars.object) = 'gwlars.object'
    return(gwlars.object)
}