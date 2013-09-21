gwselect.fit.oracle = function(x, y, coords, indx=NULL, loc, bw=NULL, family='gaussian', dist=NULL, oracle=NULL, tuning=FALSE, predict=FALSE, simulation=FALSE, verbose=FALSE, interact=FALSE, mode.select, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, N=N, AICc) {
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
	
    #Establish the oracular data frame (possibly with location interactions)
    xx = as.matrix(x[,oracle])   
    colnames(xx) = c(oracle)
    if (interact && ncol(xx)>0) {
        newnames = vector()
        for (l in 1:length(oracle)) {
            newnames = c(newnames, paste(oracle[l], ":", colnames(coords)[1], sep=""))
            newnames = c(newnames, paste(oracle[l], ":", colnames(coords)[2], sep=""))
        }

        interacted = matrix(ncol=2*ncol(xx), nrow=nrow(xx))
        for (k in 1:ncol(xx)) {
            interacted[,2*(k-1)+1] = xx[,k] * (coords[,1]-loc[1,1])
            interacted[,2*k] = xx[,k] * (coords[,2]-loc[1,2])
        }
        xx = cbind(xx, interacted)
        colnames(xx) = c(oracle, newnames)
    }

    ssr.local = NA
    df = ncol(xx) + 1
	#xx = as.matrix(x)
    yy = as.matrix(y)    
    w <- prior.weights * gwr.weights
    
    if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    xx = xx[weighted,]
    yy = as.matrix(yy[weighted])
    w = w[weighted]
    colocated = which(gwr.weights[weighted]==1)
    fitdata = data.frame(y=yy, xx)
    localdata = data.frame(fitdata[colocated,])
    colnames(fitdata) = colnames(localdata) = c("y", colnames(xx))
    
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

        model = glm(y~., data=permuted, weights=w, family=family)
        colocated = which(gwr.weights[weighted][permutation]==1)
        yyy = yy[permutation]
                
        #Get the coefficients:
        coefs = rep(0, ncol(x)+1)
        names(coefs) = c("(Intercept)", colnames(x))
        coefs[c("(Intercept)", oracle)] = coef(model)[c("(Intercept)", oracle)]
        coefs = Matrix(coefs, ncol=1)
        coef.list[[i]] = coefs
    
    	if (i==N) { 
	    	if (sum(w) > dim(permuted)[2])  {  		
				s2 = sum(w[permutation]*model$residuals^2)/(sum(w) - dim(permuted)[2] - 1)
	
				#Get standard errors of the coefficient estimates:
				se.coef = rep(0, ncol(x)+1)
                names(se.coef) = c("(Intercept)", colnames(x))
				se.coef[c(1, oracle)] = summary(model)$coefficients[,'Std. Error'][c("(Intercept)", oracle)]
				se.coef = Matrix(se.coef, ncol=1)

				#Find the local loss (for tuning bw)
				if (mode.select=='CV') {
					predictions = predict(model, newdata=localdata)
					loss.local = abs(Matrix(predictions - y[colocated], ncol=1))
				} else {
                    if (mode.select=='AIC') {penalty=2}
                    else if (mode.select=='BIC') {penalty=log(sum(w[permutation]))}

					fitted = predict(model, newdata=localdata, type='response')
					s2 = sum(w[permutation]*model$residuals**2) / (sum(w) - ncol(xx) - 1)   
				    Xh = diag(sqrt(w[permutation])) %*% as.matrix(cbind(rep(1,length(permutation)), xx))
                    H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
                    Hii = sum(H[colocated,colocated])

					if (length(colocated)>0) {
						if (!AICc) {loss.local = log(s2) + penalty*df/sum(w)}
						else {
                            loss.local = Hii
                            ssr.local = sum((w[permutation]*model$residuals**2)[colocated])
                            if (family=='gaussian') {ssr.local = sum((w[permutation]*(fitted - yy[permutation])**2)[colocated])}
                            else if (family=='poisson') {ssr.local = sum((2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - (yyy-fitted)))[colocated])}
                            else if (family=='binomial') {ssr.local = sum((2*w[permutation]*(ylogy(yyy) - yyy*log(fitted) - ylogy(1-yyy) + (1-yyy)*log(1-fitted)))[colocated])}
                        }
					} else {
						loss.local = NA
					}	
				}
			} else {
				loss.local = Inf
				fitted = NA
				s2 = NA
				se.coef = NA
				coef.list = NA
				coefs = NA
			}
        }
    }
    
    #Return the results
    if (tuning) {
        return(list(loss.local=drop(loss.local), ssr.local=ssr.local, s=NULL, sigma2=s2, nonzero=oracle, weightsum=sum(w)))
    } else if (predict) {
        return(list(loss.local=drop(loss.local), coef=coefs))
    } else if (simulation) {
        return(list(loss.local=drop(loss.local), coef=coefs, coeflist=coef.list, bw=bw, sigma2=s2, se.coef=se.coef, fitted=fitted[1], nonzero=oracle, weightsum=sum(w), s=NULL))
    } else {
        return(list(model=model, coef=coefs, coeflist=coef.list, loc=loc, bw=bw, loss.local=drop(loss.local), sigma2=s2, nonzero=oracle, sum.weights=sum(w), N=N))
    }
}
