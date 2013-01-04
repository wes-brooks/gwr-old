gwlars.fit.oracle = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, oracle=NULL, tuning=FALSE, predict=FALSE, simulation=FALSE, verbose=FALSE, interact=FALSE, mode.select, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, N=N) {
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

    vars = colnames(x)
    which.oracle = vector()
    for (v in oracle) {
        if (v %in% vars)
            which.oracle = c(which.oracle, which(vars==v))
    }
    df = length(oracle) + 1
    #x = x[,oracle]


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
        
        wo = vector()
        for (l in 1:length(which.oracle)) {
            wo = c(wo, 2*(which.oracle-1) + 1)
            wo = c(wo, 2*which.oracle)
        }
        which.oracle = wo
        oracle = colnames(x)[which.oracle]
    }

    

    xx = as.matrix(x)
    yy = as.matrix(y)    
    w <- prior.weights * gwr.weights
    
    if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    xx = xx[weighted,oracle]
    yy = as.matrix(yy[weighted])
    w = w[weighted]
    colocated = which(gwr.weights[weighted]==1)
    fitdata = data.frame(y=yy, xx)    
    localdata = data.frame(fitdata[colocated,])
    colnames(fitdata) = colnames(localdata) = c("y", oracle)
    
    int.list = list()
    coef.list = list()
    
    for (i in 1:N) {
        #Final permutation is the original ordering of the data:
        if (i==N) {            
            permutation = 1:n.weighted
        } else {
            permutation = sample(1:n.weighted, replace=TRUE)
        }

        permuted = data.frame(fitdata[permutation,])
        colnames(permuted) = colnames(fitdata)
        model = lm(y~., data=permuted, weights=w)
                
        #Get the coefficients:
        coefs = rep(0, ncol(x)+1)
        coefs[c(1, which.oracle+1)] = coef(model)
        coefs = Matrix(coefs, ncol=1)
        rownames(coefs) = c("(Intercept)", colnames(x))

        if (interact) {
            locmat = t(as.matrix(loc))
            cc = Matrix(0, nrow=(length(coefs)-1)/2, ncol=2)
            cc[,1] = coefs[seq(2, length(coefs)-1, by=2)]
            cc[,2] = coefs[seq(2, length(coefs)-1, by=2)+1]         
            ccc = cc %*% locmat
            coefs = Matrix(c(coefs[1], as.vector(ccc)), ncol=1)
            rownames(coefs) =  c("(Intercept)", oldnames)
        }


        coef.list[[i]] = coefs
        if (verbose) {print(coefs)}
    
        if (i==N) { 
            s2 = sum(model$residuals^2)/sum(w[permutation])
    
            #Get standard errors of the coefficient estimates:
            se.coef = rep(0, ncol(x)+1)
            se.coef[c(1, which.oracle+1)] = summary(model)$coefficients[,'Std. Error']
            se.coef = Matrix(se.coef, ncol=1)
            rownames(se.coef) = c("(Intercept)", colnames(x))
            
            if (interact) {
                locmat = t(as.matrix(loc**2))
                cc = Matrix(0, nrow=(length(se.coef)-1)/2, ncol=2)
                cc[,1] = se.coef[seq(2, length(se.coef)-1, by=2)]**2
                cc[,2] = se.coef[seq(2, length(se.coef)-1, by=2)+1]**2           
                ccc = sqrt(cc %*% locmat)
                se.coef = Matrix(c(se.coef[1], as.vector(ccc)), ncol=1)
                rownames(se.coef) =  c("(Intercept)", oldnames)
            } 
            
    
            #Find the local loss (for tuning bw)
            if (mode.select=='CV') {
                predictions = predict(model, newdata=localdata)
                loss.local = abs(Matrix(predictions - y[colocated], ncol=1))      
    
            } else if (mode.select=='AIC') {                           
                fitted = predict(model, newdata=localdata)
                s2 = sum(w[permutation]*model$residuals**2) / sum(w[permutation])    
                
                if (length(colocated)>0) {
                    loss.local = log(s2) + 2*df/sum(w[permutation])
                } else {
                    loss.local = NA
                }                     
    
            } else if (mode.select=='BIC') {   
                fitted = predict(model, newdata=localdata)
                s2 = sum(w[permutation]*model$residuals**2) / sum(w[permutation])
    
                if (length(colocated)>0) {
                    loss.local = log(s2) + 2*log(sum(w[permutation]))/sum(w[permutation])
                } else {
                    loss.local = NA
                }
            }
        }
    }
    
    #Return the results
    if (tuning) {
        return(list(loss.local=loss.local))
    } else if (predict) {
        return(list(loss.local=loss.local, coef=coefs))
    } else if (simulation) {
        return(list(loss.local=loss.local, coef=coefs, coeflist=coef.list, bw=bw, sigma2=s2, se.coef=se.coef))
    } else {
        return(list(model=model, coef=coefs, coeflist=coef.list, loc=loc, bw=bw, loss.local=loss.local, sigma2=s2, sum.weights=sum(w), N=N))
    }
}