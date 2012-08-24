gwglmnet.fit <-
function(x, y, coords, weight.matrix, s, verbose, family, prior.weights) {
#Fit the gwglmnet model (non-adaptive algorithm)
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    gwglmnet.object = list()
    cv.error = list()

    for(i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[,1]==coords.unique[i,1] & coords[,2]==coords.unique[i,2])
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
