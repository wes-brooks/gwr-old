gwlars.fit.inner = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, s=NULL, mode.select='', tuning=FALSE, predict=FALSE, simulation=FALSE, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, interact, precondition, N=1, tau=3, shrunk.fit) {
    colocated = which(round(coords[,1],5) == round(as.numeric(loc[1]),5) & round(coords[,2],5) == round(as.numeric(loc[2]),5))        
    reps = length(colocated)

    if (is.null(gwr.weights)) {
        gwr.weights = gweight(dist, bw)     
    } 
    gwr.weights = drop(gwr.weights) 

    #Get names for the dataset with interactions on the location variables.
    newnames = vector()
    oldnames = colnames(x)
    for (l in 1:length(oldnames)) {
        newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[1], sep=""))
        newnames = c(newnames, paste(oldnames[l], ":", colnames(coords)[2], sep=""))
    }

    #Create the dataset with interactions on the location variables.
    interacted = matrix(ncol=2*ncol(x), nrow=nrow(x))
    for (k in 1:ncol(x)) {
        interacted[,2*(k-1)+1] = x[,k]*coords[,1]
        interacted[,2*k] = x[,k]*coords[,2]
    }
    x.interacted = cbind(x, interacted)
    colnames(x.interacted) = c(oldnames, newnames)

    #Put the data in matrices
    yy = as.matrix(y)
    xx = as.matrix(x)
    xx.interacted = as.matrix(x.interacted)
    
    #Return right now if the bandwidth is so tight that there are no other observations within the kernel.
    if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

    #Find which and how many observations have nonzero weight
    n <- nrow(xx)
    w <- prior.weights * gwr.weights
    weighted = which(w>0)
    colocated = which(gwr.weights[weighted]==1)
    n.weighted = length(weighted)

    #Limit our focus to observations with nonzero weight
    w = w[weighted]
    xx = xx[weighted,]
    yy = yy[weighted]
    xx.interacted = xx.interacted[weighted,]

    #Generate a tuning dataset (with the colocated points removed)
    xx.tuning = xx[-colocated,]
    yy.tuning = yy[-colocated,]

    #Multiply by the square root of weight
    sqrt.w.tuning <- diag(sqrt(w[-colocated])) 
    sqrt.w = diag(sqrt(w)) 
    
    #Response multiplied by sqrt(weight)
    yyy = sqrt.w %*% yy
    yyy.tuning = sqrt.w.tuning %*% yy.tuning
    meany = sum((w*yy))/sum(w)
    meany.tuning = sum((w*yy)[-colocated])/sum(w[-colocated])
    
    #Covariates multiplied by sqrt(weight)
    xxx = sqrt.w %*% xx
    xxx.tuning = sqrt.w.tuning %*% xx.tuning
    xxx.interacted = sqrt.w %*% xx.interacted
    xxx.tuning.interacted = sqrt.w.tuning %*% xx.tuning.interacted

    #Add an intercept column before sending them off to lars
    tunex = cbind(sqrt(w[-colocated]), xxx.tuning)
    tuney = yyy.tuning
    fitx = cbind(sqrt(w), xxx)
    fity = yyy

    model = lars(x=tunex, y=tuney, type='lar', normalize=FALSE, intercept=FALSE)
    nsteps = length(model$lambda) + 1

    if (mode.select=='CV') {
        reps = length(colocated)
        vars = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {which(abs(x)>0)})
        df = sapply(vars, length) + 1                        

        predictions = predict(model, newx=xx[colocated,], type='fit', mode='step')[['fit']]
        loss = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=nsteps), nrow=reps, ncol=nsteps)))                
        s2 = sum(w*lsfit(y=predy, x=predx, wt=w)$residuals**2) / sum(w) 

        loss.local = loss  
        k = which.min(loss)     

    } else if (mode.select=='AIC') {
        #predx = cbind(1, t(apply(xx[permutation,], 1, function(X) {(X-meanx) * adapt.weight / normx})))
        predx = cbind(1, as.matrix(xx[permutation,]))
        predy = as.matrix(yy[permutation])

        if (sum(w) > ncol(x)) {               
            #Get the un-penalized intercept
            wcm = apply(predx[,-1], 2, function(x) {sum(x*w[permutation])})/sum(w)
            coefs = as.matrix(coef(model))
            coefs = t(apply(coefs, 1, function(x) {x * c(1, adapt.weight) * c(1, 1/normx)}))
            coefs[,1] = apply(coefs, 1, function(x) {x[1] - sum(x[-1] * wcm)})
            fitted = predx %*% t(coefs)
            coefs[,1] = coefs[,1] + apply(fitted, 2, function(x) {sum(w[permutation]*(predy-x)) / sum(w)})
            fitted = predx %*% t(coefs)    
            vars = apply(coefs, 1, function(x) {which(abs(x[-1])>0)})
            df = sapply(vars, length) + 1
            
            s2 = sum(w[permutation]*(fitted[,ncol(fitted)] - yy[permutation])**2) / (sum(w) - ncol(x))#ncol(x))  
            loss = as.vector(apply(fitted, 2, function(z) {sum(w[permutation]*(z - yy[permutation])**2)})/s2 + log(s2) + 2*df)
            k = which.min(loss)
            fitted = fitted[,k]
            localfit = fitted[colocated]      
            df = df[k]          

            if (k > 1) {
                varset = vars[[k]] #- 1
                modeldata = data.frame(y=yy[permutation], xx[permutation,varset])
                m = lm(y~., data=modeldata, weights=w[permutation])
                if (!shrunk.fit) {
                    fitted = m$fitted
                    localfit = fitted[colocated]
                    df = length(varset) + 1
                    s2 = sum((m$residuals*w[permutation])**2) / (sum(w) - df)     
                }

                coefs.unshrunk = rep(0, ncol(xx) + 1)
                coefs.unshrunk[c(1, varset + 1)] = coef(m)
                s2.unshrunk = sum(m$residuals**2)/sum(w[permutation])

                se.unshrunk = rep(0, ncol(x) + 1)
                se.unshrunk[c(1, varset + 1)] = summary(m)$coefficients[,'Std. Error']
                
                if (interact) {
                    varset.interacted = vars[[k]] #- 1
                    for (j in vars[[k]]) {
                        varset.interacted = c(varset.interacted, ncol(x)+2*(j-1)+1, ncol(x)+2*j)
                    }			

                    modeldata = data.frame(y=yy[permutation], xx.interacted[permutation,varset.interacted])
                    m = lm(y~., data=modeldata, weights=w[permutation])                        
                    if (!shrunk.fit) {
                        fitted = m$fitted
                        localfit = fitted[colocated]
                        df = length(varset.interacted) + 1
                        s2 = sum((m$residuals*w[permutation])**2) / sum(w)                            
                    }

                    coefs.unshrunk.interacted = rep(0, ncol(xx.interacted) + 1)
                    coefs.unshrunk.interacted[c(1, varset.interacted + 1)] = coef(m)
                    s2.unshrunk.interacted = sum(m$residuals**2)/sum(w[permutation])

                    se.unshrunk.interacted = rep(0, ncol(xx.interacted) + 1)
                    se.unshrunk.interacted[c(1, varset.interacted + 1)] = summary(m)$coefficients[,'Std. Error']
                }
                else {
                    coefs.unshrunk.interacted = c(meany, rep(0, ncol(xx.interacted)))
                    se.unshrunk.interacted = c(se.unshrunk[1], rep(0, ncol(xx.interacted)))
                    s2.unshrunk.interacted = s2.unshrunk
                }
            } else {
                coefs.unshrunk = rep(0, ncol(xx) + 1)
                coefs.unshrunk[1] = meany
                
                s2.unshrunk = sum(fity**2)/sum(w[permutation])
                se.unshrunk = rep(0, ncol(xx) + 1)
                se.unshrunk[1] = sqrt(s2.unshrunk)
                
                coefs.unshrunk.interacted = c(meany, rep(0, ncol(xx.interacted)))
                se.unshrunk.interacted = c(se.unshrunk[1], rep(0, ncol(xx.interacted)))
                s2.unshrunk.interacted = s2.unshrunk
            }
            
            if (length(colocated)>0) {
                loss.local = sum((w[permutation]*(fitted - yy[permutation])**2)[colocated])/s2 + log(s2) + 2*df/sum(w[permutation])
            } else {
                loss.local = NA
            }                     
        } else {
            vars = c()
            df=1
            s2 = 0
            loss = Inf
            loss.local = c(Inf)   
            fitted = rep(meany, length(permutation))
            localfit = meany
        }
    }

    #Get the tuning parameter to minimize the loss:
    s.optimal = which.min(loss)
    #loss.local = loss.local[s.optimal]

    #We have all we need for the tuning stage.
    if (!tuning) {
        #Get the coefficients:
        coefs = coefs[s.optimal,]
        coefs = Matrix(coefs, ncol=1)
        rownames(coefs) = c("(Intercept)", colnames(x))   
            
        #coefs = coefs * c(1, adapt.weight) * c(1, 1/normx)
        if (length(coefs)>1) {coefs[1] = coefs[1] - sum(coefs[2:length(coefs)] * meanx)}

        #if (interact) {
        #    locmat = t(as.matrix(cbind(rep(1,nrow(loc)),loc)))
        #    cc = Matrix(0, nrow=(length(coefs)-1-length(oldnames))/2, ncol=3)
        #    cc[,1] = coefs[seq(2, 1+length(oldnames))]
        #    cc[,2] = coefs[seq(2+length(oldnames), length(coefs)-1, by=2)]
        #    cc[,3] = coefs[seq(2+length(oldnames), length(coefs)-1, by=2)+1]            
        #    ccc = cc %*% locmat
        #    coefs = Matrix(c(coefs[1], as.vector(ccc)))
        #    rownames(coefs) =  c("(Intercept)", oldnames)
        #}     

        coefs.unshrunk = Matrix(coefs.unshrunk, ncol=1)
        rownames(coefs.unshrunk) = c("(Intercept)", colnames(xx))

        if (interact) {
            coefs.unshrunk.interacted = Matrix(coefs.unshrunk.interacted, ncol=1)
            rownames(coefs.unshrunk.interacted) = c("(Intercept)", colnames(xx.interacted))
        
            locmat = t(as.matrix(cbind(rep(1,nrow(loc)),loc)))
            cc = Matrix(0, nrow=(length(coefs.unshrunk.interacted)-1-length(oldnames))/2, ncol=3)
            cc[,1] = coefs.unshrunk.interacted[seq(2, 1+length(oldnames))]
            cc[,2] = coefs.unshrunk.interacted[seq(2+length(oldnames), length(coefs.unshrunk.interacted)-1, by=2)]
            cc[,3] = coefs.unshrunk.interacted[seq(2+length(oldnames), length(coefs.unshrunk.interacted)-1, by=2)+1]            
            ccc = cc %*% locmat
            coefs.unshrunk.interacted = Matrix(c(coefs.unshrunk.interacted[1], as.vector(ccc)))
            rownames(coefs.unshrunk.interacted) =  c("(Intercept)", oldnames)
            
            coef.unshrunk.interacted.list[[i]] = coefs.unshrunk.interacted
        }     

        se.unshrunk = Matrix(se.unshrunk, ncol=1)
        rownames(se.unshrunk) = c("(Intercept)", colnames(xx))
        
        if (interact) {
            se.unshrunk.interacted = Matrix(se.unshrunk.interacted, ncol=1)
            rownames(se.unshrunk.interacted) = c("(Intercept)", colnames(xx.interacted))
        
            locmat = t(as.matrix(cbind(rep(1,nrow(loc)),loc)))
            cc = Matrix(0, nrow=(length(se.unshrunk.interacted)-1-length(oldnames))/2, ncol=3)
            cc[,1] = se.unshrunk.interacted[seq(2, 1+length(oldnames))]**2
            cc[,2] = se.unshrunk.interacted[seq(2+length(oldnames), length(se.unshrunk.interacted)-1, by=2)]**2
            cc[,3] = se.unshrunk.interacted[seq(2+length(oldnames), length(se.unshrunk.interacted)-1, by=2)+1]**2           
            ccc = sqrt(cc %*% locmat)
            se.unshrunk.interacted = Matrix(c(se.unshrunk.interacted[1], as.vector(ccc)))
            rownames(se.unshrunk.interacted) =  c("(Intercept)", oldnames)
        }     

        coef.unshrunk.list[[i]] = coefs.unshrunk
        coef.list[[i]] = coefs
    }
    
    if (tuning) {
        return(list(loss.local=loss.local, s=s.optimal, sigma2=s2, nonzero=colnames(x)[vars[[s.optimal]]], weightsum=sum(w)))
    } else if (predict) {
        return(list(loss.local=loss.local, coef=coefs))
    } else if (simulation) {
        return(list(loss.local=loss.local, coef=coefs, coeflist=coef.list, s=s.optimal, bw=bw, sigma2=s2, coef.unshrunk=coefs.unshrunk, coef.unshrunk.interacted=coefs.unshrunk.interacted, s2.unshrunk=s2.unshrunk, s2.unshrunk.interacted=s2.unshrunk.interacted, coef.unshrunk.list=coef.unshrunk.list, coef.unshrunk.interacted.list=coef.unshrunk.interacted.list, se.unshrunk=se.unshrunk, se.unshrunk.interacted=se.unshrunk.interacted, nonzero=colnames(x)[vars[[s.optimal]]], fitted=localfit, weightsum=sum(w), loss=loss))
    } else {
        return(list(model=model, loss=loss, coef=coefs, coef.unshrunk=coefs.unshrunk, coeflist=coef.list, s=s.optimal, loc=loc, bw=bw, meanx=meanx, meany=meany, coef.scale=adapt.weight/normx, df=df, loss.local=loss.local, sigma2=s2, sum.weights=sum(w), N=N, fitted=localfit, coef.unshrunk.interacted.list=coef.unshrunk.interacted.list, s2.unshrunk.interacted=s2.unshrunk.interacted, coef.unshrunk.interacted=coefs.unshrunk.interacted, se.unshrunk.interacted=se.unshrunk.interacted))
    }
}
