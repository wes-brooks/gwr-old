gwlars.fit.inner = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, s=NULL, mode.select='', tuning=FALSE, predict=FALSE, simulation=FALSE, shrink=TRUE, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, mode, precondition=FALSE, N=N) {
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
    coef.unshrunk.list=list()

    for (i in 1:N) {
        #Final permutation is the original ordering of the data:
        if (i==N) {            
            permutation = 1:n.weighted
        } else {
            permutation = sample(1:n.weighted, replace=TRUE)
        }

        sqrt.w <- diag(sqrt(w[permutation]))       
    
        xxx = sqrt.w %*% xx[permutation,]
        yyy = sqrt.w %*% yy[permutation,]
        meany = mean(yyy)
        yyy = yyy - meany

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
            
            lm.step = try(lm(yyy~xs-1))  # mle fit on standardized
        
            if(class(lm.step) == "try-error") { 
                cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
                return(return(list(loss.local=Inf, resid=Inf)))
            }
            
            beta.lm = lm.step$coeff #[2:(m+1)]                   # mle except for intercept
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
    
        model = lars(x=xfit, y=yfit, type='lar', normalize=TRUE, intercept=FALSE)
        nsteps = length(model$lambda) + 1
    
        if (mode.select=='CV') {
            predx = matrix(predx, reps, dim(xs)[2])    
            vars = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {which(abs(x)>0)})
            df = sapply(vars, length) + 1            
            
            if (shrink==FALSE) {
                models = list()
                loss = vector()
                
                for (k in 1:nsteps) {
                    varset = vars[[k]]
                    modeldata = cbind(yfit, xfit[,varset])
                    models[[k]] = lm(yfit~., data=modeldata, weights=w)         
                    loss = c(loss, sum(abs(predict(models[[k]], newx=predx)-y[colocated])))
                    rownames(cc)=names(coef(mm))
                }
            } else {
                predictions = predict(model, newx=predx, type='fit', mode='step')[['fit']]
                loss = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=nsteps), nrow=reps, ncol=nsteps)))                
                s2 = sum(w*(fitted[,nsteps] - as.matrix(y))**2) / (sum(w)-df-1)
            }
            loss.local = loss        

        } else if (mode.select=='AIC') {
            predx = t(apply(xx, 1, function(X) {(X-meanx) * adapt.weight / normx}))
            vars = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {which(abs(x)>0)})
            df = sapply(vars, length) + 1
                        
            if (sum(w[permutation]) > nsteps) {
                if (shrink==FALSE) {
                    cv = cv.lars(x=xfit, y=yfit, mode='step', type='lar', K=min(10, n.weighted), plot.it=FALSE)
                    cv.lim = min(cv$cv) + cv$cv.error[which.min(cv$cv)]
                    S = which(cv$cv<cv.lim)[1]

                    models = list()
                    loss = vector()
                    coefs = vector()
                    fitted = vector()
                    s2 = vector()
                    for (k in 1:nsteps) {
                        varset = vars[[k]]
                        modeldata = data.frame(y=yy[permutation], predx[permutation,varset])
                        colnames(modeldata) = c('y', names(varset))
                        m = lm(y~.-1, data=modeldata, weights=w)
                        coefs = rbind(coefs, rep(0, ncol(x)+1))
                        coefs[nrow(coefs), c(1,varset+1)] = c(meany, coef(m))
                        fitted = cbind(fitted, predict(m, newx=predx))
                        loss = c(loss, ifelse(k==S,0,Inf))
                        s2 = c(s2, summary(m)$sigma**2)
                    }
                    #print(loss)
                } else {
                    coefs = cbind(`(Intercept)` = meany, coef(model))
                    fitted = predict(model, newx = predx, type = "fit", mode = "step")[["fit"]]
                    s2 = sum(lsfit(y=yfit, x=xfit)$residuals^2) / (sum(w[permutation]) - nsteps - 1)
                    loss = as.vector(apply(fitted, 2, function(z) {sum((w * (z - yy))[permutation]^2)}) / s2 + 2 * df)
                    k = which.min(loss)
                
                    if (k > 1) {
                        varset = vars[[k]]
                        modeldata = data.frame(y=yy[permutation], xx[permutation,varset])
                        m = lm(y~., data=modeldata, weights=w)
                        #m = lsfit(x=modeldata, y=yfit, intercept=FALSE)
                        coefs.unshrunk = rep(0, ncol(x) + 1)
                        coefs.unshrunk[c(1, varset + 1)] = coef(m)
                        s2.unshrunk = sum(m$residuals^2)/(sum(w[permutation]) - 1 - length(coef(m)))

                        se.unshrunk = rep(0, ncol(x) + 1)
                        se.unshrunk[c(1, varset + 1)] = summary(m)$coefficients[,'Std. Error']
                    } else {
                        coefs.unshrunk = rep(0, ncol(x) + 1)
                        coefs.unshrunk[1] = meany
                        
                        s2.unshrunk = sum(yfit^2)/(sum(w[permutation]) - 1)
                        se.unshrunk = rep(0, ncol(x) + 1)
                        se.unshrunk[1] = sqrt(s2.unshrunk / sum(w[permutation]))
                    }
                }

                if (length(colocated)>0) {
                    loss.local = as.vector(apply(fitted, 2, function(z) {sum((w*(z - (yy-meany)))[colocated]**2)})/s2 + log(s2) + 2*df/sum(w[permutation])) 
                } else {
                    loss.local = rep(NA, length(loss))
                }                     
                #print(loss.local)
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
                loss.local = as.vector(apply(fitted, 2, function(z) {sum((w*(z - y))[colocated]**2)})/s2 + log(s2) + log(sum(w[permutation]))*df/sum(w[permutation]))
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
        coefs = coefs[s.optimal,] #predict(model, type='coefficients', s=s.optimal, mode='step')[['coefficients']]
        coefs = Matrix(coefs, ncol=1)
        rownames(coefs) = c("(Intercept)", colnames(x))
    
        coefs = coefs * c(1, adapt.weight) / c(1, normx)
        coefs[1] = coefs[1] - sum(coefs[2:length(coefs)] * meanx)
        #print(coefs)
        if (verbose) {print(coefs)}

        coefs.unshrunk = Matrix(coefs.unshrunk, ncol = 1)
        rownames(coefs.unshrunk) = c("(Intercept)", colnames(x))

        se.unshrunk = Matrix(se.unshrunk, ncol = 1)
        rownames(se.unshrunk) = c("(Intercept)", colnames(x))
        
        #coefs.unshrunk = coefs.unshrunk * c(1, adapt.weight) / c(1, normx)
        #coefs.unshrunk[1] = coefs.unshrunk[1] - sum(coefs.unshrunk[2:length(coefs.unshrunk)] * meanx)

        coef.unshrunk.list[[i]] = coefs.unshrunk
        coef.list[[i]] = coefs
    }
    
        
    #Get the residuals at this choice of s:
    #fitted = predict(model, newx=xfit, s=s.optimal, type='fit', mode='step')[['fit']]
    #resid = yfit - fitted
    
    if (tuning) {
        return(list(loss.local=loss.local))
    } else if (predict) {
        return(list(loss.local=loss.local, coef=coefs))
    } else if (simulation) {
        return(list(loss.local=loss.local, coef=coefs, coeflist=coef.list, s=s.optimal, bw=bw, sigma2=s2, coefs.unshrunk=coefs.unshrunk, s2.unshrunk=s2.unshrunk, coef.unshrunk.list=coef.unshrunk.list, se.unshrunk=se.unshrunk))
    } else {
        return(list(model=model, loss=loss, coef=coefs, coeflist=coef.list, s=s.optimal, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, df=df, loss.local=loss.local, sigma2=s2, sum.weights=sum(w), N=N))
    }
}
