gwlars.fit.inner = function(x, y, coords, loc, bw=NULL, dist=NULL, s=NULL, mode.select='', verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, mode, precondition=FALSE, N=N) {
    colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))
    reps = length(colocated)

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } else {
        gwr.weights = gwr.weights
    }      

    if (sum(gwr.weights)==0) { return(list(cv.error=Inf, resid=Inf)) }   
      
    if (mode.select=='CV') { 
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
        w <- prior.weights[-colocated] * gwr.weights[-colocated]   
    } else {
        xx = as.matrix(x)
        yy = as.matrix(y)
        w <- prior.weights * gwr.weights
    }

    m <- ncol(xx)
    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    int.list = list()
    coef.list = list()

    for (i in 1:N) {
        #Final permutation is the original ordering of the data:
        if (i==N) {            
            permutation = 1:n    
        } else {
            permutation = sample(weighted, replace=TRUE)
        }

        sqrt.w <- diag(sqrt(w[permutation]))       
    
        xxx = sqrt.w %*% xx[permutation,]
        yyy = sqrt.w %*% yy[permutation,]
        
        if (precondition==TRUE) {
            s = svd(xxx)
            F = s$u  %*% diag(1/s$d)  %*%  t(s$u)
            xxx = F %*% xxx
            yyy = F %*% yyy
        }
    
        if (adapt==TRUE) {
            one <- rep(1, nrow(xxx))
            meanx <- drop(one %*% xxx)/nrow(xxx)
            x.centered <- scale(xxx, meanx, FALSE)         # first subtracts mean
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
            
            lm.step = try(lm(yyy~xs, weights=w[permutation]))  # mle fit on standardized
        
            if(class(lm.step) == "try-error") { 
                cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
                return(Inf)
            }
            
            beta.lm = lm.step$coeff[2:(m+1)]                   # mle except for intercept
            adapt.weight = abs(beta.lm)                        # weights for adaptive lasso
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(adapt.weight[k])) {
                    xs[,k] = xs[,k] * adapt.weight[k]
                } else {
                    xs[,k] = rep(0, dim(xs)[1])
                    adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
                }
            }
            predx = (x[colocated,] - meanx) * adapt.weight / normx
            
        } else {
            meanx = rep(0, ncol(x))
            adapt.weight = rep(1, ncol(x))
            normx = rep(1, ncol(x))
    
            xs=xxx
            predx = x[colocated,]
        }
    
        xfit = xs #sqrt.w %*% xs
        yfit = yyy #sqrt.w %*% yyy
    
        model = lars(x=xfit, y=yfit, type='lar', normalize=TRUE, intercept=TRUE)
        #ll = model$lambda
        nsteps = length(model$lambda) + 1
    
        if (mode.select=='CV') {
            df=NULL
            predx = matrix(predx, reps, dim(xs)[2])
            #predictions = predict(model, newx=predx, s=ll, type='fit', mode='step')[['fit']]
            predictions = predict(model, newx=predx, type='fit', mode='step')[['fit']]
            #cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(ll)), nrow=reps, ncol=length(ll))))
            loss = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=nsteps), nrow=reps, ncol=nsteps)))
            loss.local = loss
    
            df = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {sum(abs(x)>0)}) + 1
            s2 = sum(w*(fitted[,nsteps] - as.matrix(y))**2) / (sum(w)-df-1)
            #s.optimal = ll[which.min(cv.error)]
        } else if (mode.select=='AIC') {
            predx = t(apply(x, 1, function(X) {(X-meanx) * adapt.weight / normx}))
            coef = predict(model, type='coefficients')[['coefficients']]
            df = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {sum(abs(x)>0)}) + 1
            fitted = predict(model, newx=predx, type='fit', mode='step')[['fit']]

            if (sum(w[permutation]) > nsteps) {
                s2 = sum((w*(fitted[,nsteps] - as.matrix(y)))[permutation]**2) / sum(w[permutation])
                loss = as.vector(apply(fitted, 2, function(z) {sum((w*(z - y))[permutation]**2)})/s2 + 2*df) 
                loss.local = as.vector(apply(fitted, 2, function(z) {sum((w*(z - y))[colocated]**2)})/s2 + log(s2) + 2*df/sum(w[permutation]))
            } else {
                s2 = 0
                loss = Inf
                loss.local = c(Inf)
            }
                
        } else if (mode.select=='BIC') {
            predx = t(apply(x[permutation,], 1, function(X) {(X-meanx) * adapt.weight / normx}))
            coef = predict(model, type='coefficients')[['coefficients']]
            df = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {sum(abs(x)>0)}) + 1
            fitted = predict(model, newx=predx, type='fit', mode='step')[['fit']]

            if (sum(w) > nsteps) {
                s2 = sum(w*(fitted[,nsteps] - as.matrix(yyy))**2) / sum(w)
                loss = as.vector(apply(fitted, 2, function(z) {sum(w*(z - yyy)**2)})/s2 + log(sum(w))*df) 
                loss.local = as.vector(apply(fitted, 2, function(z) {sum((w*(z - y))[colocated]**2)})/s2 + log(s2) + log(sum(w))*df/sum(w))
            } else {
                s2 = 0
                loss = Inf
                loss.local = c(Inf)
            }
        }
    
        #Get the tuning parameter to minimize the loss:
        s.optimal = which.min(loss)
        loss.local = loss.local[s.optimal]
        
        #Get the coefficients:
        coef = predict(model, type='coefficients', s=s.optimal, mode='step')[['coefficients']]
        coef = Matrix(coef, ncol=1)
        rownames(coef) = colnames(x)
    
        coef = coef * adapt.weight / normx
        intercept = predict(model, type='fit', s=s.optimal, mode='step', newx=matrix(0,1,nrow(coef)))[['fit']] #- coef[2:length(coef)] * meanx[2:length(coef)]

        int.list[[i]] = intercept
        coef.list[[i]] = coef
    }
        
    #Get the residuals at this choice of s:
    fitted = predict(model, newx=xfit, s=s.optimal, type='fit', mode='step')[['fit']]
    resid = yfit - fitted
    
    return(list(model=model, loss=loss, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=resid, coef=coef, intercept=intercept, coeflist=coef.list, intlist=int.list, df=df, loss.local=loss.local, sigma2=s2, sum.weights=sum(w), N=N))
}