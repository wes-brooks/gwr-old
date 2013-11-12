gwglmnet.fit.ridge = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, s=NULL, family, verbose, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, interact, precondition, tau=3) {
    if (!is.null(indx)) {
        colocated = which(round(coords[indx,1],5)==round(as.numeric(loc[1]),5) & round(coords[indx,2],5) == round(as.numeric(loc[2]),5))
    }
    else {
        colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))        
    }
    reps = length(colocated)

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } 
    gwr.weights = drop(gwr.weights)  

    if (!is.null(indx)) {
        gwr.weights = gwr.weights[indx]
    }

    #For interaction on location:
    if (interact) {
        newnames = vector()
        oldnames = colnames(x)
        for (l in 1:length(oldnames)) {
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[1], sep=""))
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[2], sep=""))
        }
        interacted = matrix(ncol=2*ncol(x), nrow=nrow(x))
        for (k in 1:ncol(x)) {
            interacted[,2*(k-1)+1] = x[,k]*(coords[,1]-loc[1,1])
            interacted[,2*k] = x[,k]*(coords[,2]-loc[1,2])
        }
        x.interacted = cbind(x, interacted)
        colnames(x.interacted) = c(oldnames, newnames)
    }


    xx = as.matrix(x)
    if (interact) {xx.interacted = as.matrix(x.interacted)}
    yy = as.matrix(y)
    w <- prior.weights * gwr.weights

    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    xx = xx[weighted,]
    if (interact) {xx.interacted = xx.interacted[weighted,]}
    yy = as.matrix(yy[weighted])
    w = w[weighted]
    
    colocated = which(gwr.weights[weighted]==1)
    sqrt.w <- diag(sqrt(w))

    if (interact) { 
        xxx = xx.interacted
    } else {
        xxx = xx
    }
    
    yyy = yy
    meany = sum(w*yy)/sum(w)

    if (precondition==TRUE) {
        s = svd(xxx)
        F = s$u  %*% diag(1/sqrt(s$d**2 + tau))  %*%  t(s$u)
        xxx = F %*% xxx
        yyy = F %*% yyy
    }
    
    if (family == 'binomial') { ridge = glmnet(x=xxx, y=cbind(1-yyy, yyy), standardize=FALSE, intercept=TRUE, family=family, weights=w, alpha=0) }
    else { ridge = glmnet(x=xxx, y=yyy, standardize=FALSE, intercept=TRUE, family=family, weights=w, alpha=0) }
    nsteps = length(ridge$lambda) + 1
        
        
}
