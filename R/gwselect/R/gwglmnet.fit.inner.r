gwglmnet.fit.inner = function(x, y, family, coords, loc, bw=NULL, tuning=FALSE, predict=FALSE, dist=NULL, indx=NULL, s=NULL, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, mode, mode.select, gweight=NULL, shrink=TRUE, longlat=FALSE, adapt=FALSE, precondition=FALSE, N) {
    if (!is.null(indx)) {
        colocated = which(coords[indx,1]==as.numeric(loc[1]) & coords[indx,2]==as.numeric(loc[2]))
    }
    else {
        colocated = which(coords[,1]==as.numeric(loc[1]) & coords[,2]==as.numeric(loc[2]))
    }
    
    reps = length(colocated)   

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } else {
        gwr.weights = gwr.weights
    }      
    gwr.weights = drop(gwr.weights)

    if (!is.null(indx)) {
        gwr.weights = gwr.weights[indx]
    }
    
    if (mode.select=='CV') { 
        xx = as.matrix(x[-colocated,])
        yy = as.matrix(y[-colocated])
        w <- prior.weights[-colocated] * gwr.weights[-colocated]  
    } else {
        xx = as.matrix(x)
        yy = as.matrix(y)
        w <- prior.weights * gwr.weights
    }

    if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

    
    m <- ncol(xx)
    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)

    xx = xx[weighted,]
    yy = as.matrix(yy[weighted])

    w = w[weighted]
    colocated = which(gwr.weights[weighted]==1)

    int.list = list()
    coef.list = list()

    for (i in 1:N) {
        #Final permutation is the original ordering of the data:
        if (i==N) {            
            permutation = 1:n.weighted
        } else {
            permutation = sample(1:n.weighted, replace=TRUE)
        }      
    
        xxx = xx[permutation,] #sqrt.w %*% xx[permutation,]
        yyy = yy[permutation,] #sqrt.w %*% yy[permutation,]
        meany = mean(yyy)
        #yyy = yyy - meany

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

            glm.step = try(glm(yyy~xs, weights=w[permutation], family=family))  # mle fit on standardized
    
            if("try-error" %in% class(glm.step)) { 
                cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
                return(return(list(loss.local=Inf, resid=Inf)))
            }
    
            beta.glm = glm.step$coeff[2:(m+1)]                    # mle except for intercept
            adapt.weight = abs(beta.glm)                        # weights for adaptive lasso
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
    
        #model = glmnet(x=xfit, y=yfit, family=family, weights=w[permutation], standardize=TRUE)
        #nsteps = length(model$lambda) + 1

        if (family=='binomial' && (abs(sum(yfit*w)-sum(w))<1e-4 || sum(yfit*w)<1e-4)) {
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yfit*w), ", sum of weights=", sum(w), "\n", sep=''))
            return(list(model=NULL, cv.error=0, s=Inf, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=Inf))
        } else if (family=='binomial') {
            if (mode.select=='AIC') {
                model = glmnet(x=xfit, y=cbind(1-yfit,yfit), family=family, weights=w, lambda=s)
                predx = t(apply(xx, 1, function(X) {(X-meanx) * adapt.weight / normx}))
                vars = apply(predict(model, type='coef'), 2, function(x) {which(abs(x)>0)})
                df2 = sapply(vars, length)
                            
                if (sum(w[permutation]) > ncol(x)+1) {                    
                    coefs = t(as.matrix(coef(model)))
                    fitted = predict(model, newx=predx, type='response')
                    #s2 = sum((w*(fitted[,nsteps] - as.matrix(yy)))[permutation]**2) / sum(w[permutation])
                    #s2 = sum(lsfit(y=yfit, x=xfit)$residuals**2) / (sum(w[permutation]) - nsteps - 1)
                    loss = as.vector(apply(fitted, 2, function(z) {sum((w*(z - yy))[permutation]**2/(z*(1-z)))}) + 2*df2/sum(w[permutation]))                                   
                    
                    if (length(colocated)>0) {
                        loss.local = as.vector(apply(fitted, 2, function(z) {sum(((w*(z - yy))**2 / (z*(1-z)))[colocated])}) + 2*df2/sum(w[permutation]))
                    } else {
                        loss.local = rep(NA, length(loss))
                    }                    
                } else {
                    s2 = 0
                    loss = Inf
                    loss.local = c(Inf)   
                }
    
            }
    
        } else if (family=='poisson') {
            model = glmnet(x=xfit, y=yfit, family=family, weights=w, lambda=s)
            predictions = predict(model, newx=predx, type='response')
            cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=length(model[['lambda']])), nrow=reps, ncol=length(model[['lambda']]))))
            s.optimal = model[['lambda']][which.min(cv.error)]
            V = function(mu) {mu}
        }

        #Get the tuning parameter to minimize the loss:
        s.optimal = which.min(loss)
        loss.local = loss.local[s.optimal]

        #Get the coefficients:
        coefs = coefs[s.optimal,] #predict(model, type='coefficients', s=s.optimal, mode='step')[['coefficients']]
        coefs = Matrix(coefs, ncol=1)
        rownames(coefs) = c("(Intercept)", colnames(x))
    
        coefs = coefs * c(1, adapt.weight) / c(1, normx)
        coefs[1] = coefs[1] - sum(coefs[2:length(coefs)] * meanx)
        #print(coefs)
        if (verbose) {print(coefs)}
        #intercept = predict(model, type='fit', s=s.optimal, mode='step', newx=matrix(0,1,nrow(coef)))[['fit']] #- coef[2:length(coef)] * meanx[2:length(coef)]

        #int.list[[i]] = intercept
        coef.list[[i]] = coefs
    }
    
        
    #Get the residuals at this choice of s:
    #fitted = predict(model, newx=xfit, s=s.optimal, type='fit', mode='step')[['fit']]
    #resid = yfit - fitted
    
    if (tuning) {
        return(list(loss.local=loss.local))
    } else if (predict) {
        return(list(loss.local=loss.local, coef=coefs))
    } else {return(list(model=model, loss=loss, coef=coefs, coeflist=coef.list, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, df=df2, loss.local=loss.local, sum.weights=sum(w), N=N)) }#, resid=resid, intercept=intercept, intlist=int.list))
}
