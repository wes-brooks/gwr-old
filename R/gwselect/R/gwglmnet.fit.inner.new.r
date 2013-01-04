gwlars.fit.inner = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, s=NULL, mode.select='', tuning=FALSE, predict=FALSE, simulation=FALSE, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, interact=FALSE, precondition=FALSE, N=1) {
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

    if (interact) {
        newnames = vector()
        oldnames = colnames(x)
        for (l in 1:length(oldnames)) {
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[1], sep=""))
            newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[2], sep=""))
        }

        interacted = matrix(ncol=2*ncol(x), nrow=nrow(x))
        for (k in 1:ncol(x)) {
            interacted[,2*(k-1)+1] = x[,k]*coords[,1]
            interacted[,2*k] = x[,k]*coords[,2]
        }
        x = interacted
        colnames(x) = newnames
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

    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    xx = xx[weighted,]
    yy = as.matrix(yy[weighted])
    w = w[weighted]

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

        colocated = which(gwr.weights[weighted][permutation]==1)
        sqrt.w <- diag(sqrt(w[permutation]))        
        yyy = yy[permutation,]
        meany = mean(yyy)
        #yyy = yyy - meany   
        #normy = sqrt(sum(yyy**2))
        #yyy = yyy / normy
     
        xxx = xx[permutation,]

        if (precondition==TRUE) {
            s = svd(xxx)
            F = s$u  %*% diag(1/s$d)  %*%  t(s$u)
            xxx = F %*% xxx
            yyy = F %*% yyy
        }
    
        if (adapt==TRUE) {
            one <- rep(1, nrow(xxx))
            meanx <- drop(one %*% xxx) / nrow(xxx)
            x.centered <- scale(xxx, meanx, FALSE)         # first subtracts mean
            normx <- sqrt(drop(one %*% (x.centered**2)))
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
            
            glm.step = try(glm(yyy~xs, family=family, weights=w[permutation]))  # mle fit on standardized
        
            if("try-error" %in% class(glm.step)) { 
                cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
                return(return(list(loss.local=Inf, resid=Inf)))
            }
            
            beta.glm = glm.step$coeff[-1]                   # mle except for intercept
            adapt.weight = abs(beta.glm)               # weights for adaptive lasso
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(adapt.weight[k])) {
                    xs[,k] = xs[,k] * adapt.weight[k]
                } else {
                    xs[,k] = rep(0, dim(xs)[1])
                    adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
                }
            }
            predx = as.matrix((xx[permutation,] - meanx) * adapt.weight / normx)
            
        } else {
            meanx = rep(0, ncol(x))
            adapt.weight = rep(1, ncol(x))
            normx = rep(1, ncol(x))
    
            xs=xxx
            predx = xx[permutation,]
        }
    
        fitx = xs
        fity = yyy
    
        #model = glmnet(x=fitx, y=fity, type='lar', normalize=TRUE, intercept=FALSE)
        #nsteps = length(model$lambda) + 1   
    
    
    	if (family=='binomial' && (abs(sum(yfit*w[permutation])-sum(w[permutation]))<1e-4 || sum(yfit*w[permutation])<1e-4)) {
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yfit*w[permutation]), ", sum of weights=", sum(w[permutation]), "\n", sep=''))
            return(list(model=NULL, loss.local=0, s=Inf, loc=loc, bw=bw, meanx=meanx, coef.scale=adapt.weight/normx, resid=Inf))
        } else if (family=='binomial') {
            if (mode.select=='AIC') {
                model = glmnet(x=xfit, y=cbind(1-yfit,yfit), family=family, weights=w[permutation], lambda=s)
                predx = t(apply(xx, 1, function(X) {(X-meanx) * adapt.weight / normx}))
                vars = apply(predict(model, type='coef'), 2, function(x) {which(abs(x)>0)})
                df2 = sapply(vars, length)
                            
                if (sum(w[permutation]) > ncol(x)+1) {                    
                    coefs = t(as.matrix(coef(model)))
                    fitted = predict(model, newx=predx, type='response')
                    #s2 = sum((w*(fitted[,nsteps] - as.matrix(yy)))[permutation]**2) / sum(w[permutation])
                    #s2 = sum(lsfit(y=yfit, x=xfit)$residuals**2) / (sum(w[permutation]) - nsteps - 1)
                    loss = as.vector(apply(fitted, 2, function(z) {-2*sum(w[permutation]*(yy*log(z) + (1-yy)*log(1-z))}) + 2*df2)
                    
                    if (length(colocated)>0) {
                        loss.local = as.vector(apply(fitted, 2, function(z) {-2*sum((w[permutation]*(yy*log(z) + (1-yy)*log(1-z))[colocated])}) + 2*df2/sum(w[permutation]))
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
    
    
        if (mode.select=='CV') {
            predx = matrix(xx[colocated,], reps, dim(xs)[2])    
            vars = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {which(abs(x)>0)})
            df = sapply(vars, length) + 1                        

            predictions = predict(model, newx=predx, type='fit', mode='step')[['fit']]
            loss = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=nsteps), nrow=reps, ncol=nsteps)))                
            s2 = sum(w*(fitted[,nsteps] - as.matrix(y))**2) / (sum(w)-df-1)

            loss.local = loss        

        } else if (mode.select=='AIC') {
        	
            predx = t(apply(xx[permutation,], 1, function(X) {(X-meanx) * adapt.weight / normx}))
            predy = as.matrix(yy[permutation])
            vars = apply(as.matrix(predict(model, type='coef')[['coefficients']]), 1, function(x) {which(abs(x)>0)})
            df = sapply(vars, length) + 1

            if (n.weighted > ncol(x)) {               
                coefs = cbind("(Intercept)"=meany, coef(model))
                fitted = predict(model, newx=predx, type="fit", mode="step")[["fit"]]

                #Mimics how lm() finds the error variance for weighted least squares, but uses sum(w) instead of n for the number of observations.
                #mod = lm(fity~fitx)                
                #r = mod$residuals / sqrt(w[permutation])
                #f = mod$fitted
                #wm = sum(w[permutation] * f)/sum(w)
                #mss = sum(w[permutation] * (f-wm)**2)
                #rss = sum(w[permutation] * r**2)                
                #s2 = rss / sum(w[permutation])    
                s2 = sum(w[permutation]*lsfit(y=predy, x=predx, wt=w[permutation])$residuals**2) / sum(w[permutation])  
                #s2 = summary(lm(yy[permutation]~xx[permutation,], weights=w[permutation]))$sigma**2

                #use wlars:
                #m2 = wlars(x=x, y=y, W=W)
                #sm2 = summary(m2)
                #sm2$Rss / rev(sm2$Rss)[1] + 2*sm2$Df

                loss = as.vector(apply(fitted, 2, function(z) {sum(w[permutation]*(z - yy[permutation])**2)})/s2 + 2*df)
                #loss = as.vector(apply(fitted, 2, function(z) {normy**2 * sum(w[permutation]*(z - fity)**2)})/s2 + 2*df)
                k = which.min(loss)

                if (k > 1) {
                    varset = vars[[k]]
                    modeldata = data.frame(y=yy[permutation], xx[permutation,varset])
                    m = lm(y~., data=modeldata, weights=w)
                    coefs.unshrunk = rep(0, ncol(x) + 1)
                    coefs.unshrunk[c(1, varset + 1)] = coef(m)
                    s2.unshrunk = sum(m$residuals^2)/(sum(w[permutation]) - 1 - length(coef(m)))

                    se.unshrunk = rep(0, ncol(x) + 1)
                    se.unshrunk[c(1, varset + 1)] = summary(m)$coefficients[,'Std. Error']
                } else {
                    coefs.unshrunk = rep(0, ncol(xx) + 1)
                    coefs.unshrunk[1] = meany
                    
                    s2.unshrunk = sum(fity**2)/(sum(w[permutation]) - 1)
                    se.unshrunk = rep(0, ncol(xx) + 1)
                    se.unshrunk[1] = sqrt(s2.unshrunk / sum(w[permutation]))
                }
                
                if (length(colocated)>0) {
                    loss.local = as.vector(apply(fitted, 2, function(z) {sum((w[permutation]*(z - yy[permutation])**2)[colocated])})/s2 + log(s2) + 2*df/sum(w[permutation]))
                    #loss.local = as.vector(apply(fitted, 2, function(z) {normy**2 * sum((w*(z - fity)**2)[colocated])})/s2 + 2*log(normy) + log(s2) + 2*df/sum(w[permutation])) 
                } else {
                    loss.local = rep(NA, length(loss))
                }                     
            } else {
                s2 = 0
                loss = Inf
                loss.local = c(Inf)   
            }

        } else if (mode.select=='BIC') {
            #predx = t(apply(xx, 1, function(X) {(X-meanx) * adapt.weight / normx}))
            vars = apply(as.matrix(predict(model, type='coef')[['coefficients']]), 1, function(x) {which(abs(x)>0)})
            df = sapply(vars, length) + 1

            if (sum(w[permutation]) > nsteps) {               
                coefs = cbind("(Intercept)"=meany, coef(model))
                fitted = predict(model, newx=fitx, type="fit", mode="step")[["fit"]]
                s2 = sum(lsfit(y=fity, x=fitx)$residuals**2) / (sum(w[permutation]) - nsteps - 1)     
                loss = as.vector(apply(fitted, 2, function(z) {sum(((z - fity)**2)[permutation])})/s2 + log(sum(w[permutation]))*df)
                k = which.min(loss)

                if (k > 1) {
                    varset = vars[[k]]
                    modeldata = data.frame(y=yy[permutation], xx[permutation,varset])
                    m = lm(y~., data=modeldata, weights=w)
                    coefs.unshrunk = rep(0, ncol(x) + 1)
                    coefs.unshrunk[c(1, varset + 1)] = coef(m)
                    s2.unshrunk = sum(m$residuals^2)/(sum(w[permutation]) - 1 - length(coef(m)))

                    se.unshrunk = rep(0, ncol(x) + 1)
                    se.unshrunk[c(1, varset + 1)] = summary(m)$coefficients[,'Std. Error']
                } else {
                    coefs.unshrunk = rep(0, ncol(xx) + 1)
                    coefs.unshrunk[1] = meany
                    
                    s2.unshrunk = sum(fity**2)/(sum(w[permutation]) - 1)
                    se.unshrunk = rep(0, ncol(xx) + 1)
                    se.unshrunk[1] = sqrt(s2.unshrunk / sum(w[permutation]))
                }
                
                if (length(colocated)>0) {
                    loss.local = as.vector(apply(fitted, 2, function(z) {sum(((z - fity)**2)[colocated])})/s2 + log(s2) + log(sum(w[permutation]))*df/sum(w[permutation])) 
                    #loss.local = as.vector(apply(fitted, 2, function(z) {sum((w*(z - (yy-meany)/normy)**2)[colocated])})/s2 + log(s2) + 2*df/sum(w[permutation]))
                } else {
                    loss.local = rep(NA, length(loss))
                }                     
            } else {
                s2 = 0
                loss = Inf
                loss.local = c(Inf)   
            }
        }

    
        #Get the tuning parameter to minimize the loss:
        s.optimal = which.min(loss)
        loss.local = loss.local[s.optimal]

        #We have all we need for the tuning stage.
        if (!tuning) {
            #Get the coefficients:
            coefs = coefs[s.optimal,]
            coefs = Matrix(coefs, ncol=1)
            rownames(coefs) = c("(Intercept)", colnames(x))   
                
            coefs = coefs * c(1, adapt.weight) * c(1, 1/normx)
            if (length(coefs)>1) {coefs[1] = coefs[1] - sum(coefs[2:length(coefs)] * meanx)}
            if (verbose) {print(coefs)}
    
            if (interact) {
                locmat = t(as.matrix(loc))
                cc = Matrix(0, nrow=(length(coefs)-1)/2, ncol=2)
                cc[,1] = coefs[seq(2, length(coefs)-1, by=2)]
                cc[,2] = coefs[seq(2, length(coefs)-1, by=2)+1]            
                ccc = cc %*% locmat
                coefs = Matrix(c(coefs[1], as.vector(ccc)))
                rownames(coefs) =  c("(Intercept)", oldnames)
            }     
    
            coefs.unshrunk = Matrix(coefs.unshrunk, ncol=1)
            rownames(coefs.unshrunk) = c("(Intercept)", colnames(xx))
    
            if (interact) {
                locmat = t(as.matrix(loc))
                cc = Matrix(0, nrow=(length(coefs.unshrunk)-1)/2, ncol=2)
                cc[,1] = coefs.unshrunk[seq(2, length(coefs.unshrunk)-1, by=2)]
                cc[,2] = coefs.unshrunk[seq(2, length(coefs.unshrunk)-1, by=2)+1]            
                ccc = cc %*% locmat
                coefs.unshrunk = Matrix(c(coefs.unshrunk[1], as.vector(ccc)))
                rownames(coefs.unshrunk) =  c("(Intercept)", oldnames)
            }     
    
            se.unshrunk = Matrix(se.unshrunk, ncol=1)
            rownames(se.unshrunk) = c("(Intercept)", colnames(xx))
            
            if (interact) {
                locmat = t(as.matrix(loc**2))
                cc = Matrix(0, nrow=(length(se.unshrunk)-1)/2, ncol=2)
                cc[,1] = se.unshrunk[seq(2, length(se.unshrunk)-1, by=2)]**2
                cc[,2] = se.unshrunk[seq(2, length(se.unshrunk)-1, by=2)+1]**2           
                ccc = sqrt(cc %*% locmat)
                se.unshrunk = Matrix(c(se.unshrunk[1], as.vector(ccc)))
                rownames(se.unshrunk) =  c("(Intercept)", oldnames)
            }     
    
            coef.unshrunk.list[[i]] = coefs.unshrunk
            coef.list[[i]] = coefs
        }
    }
    
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
