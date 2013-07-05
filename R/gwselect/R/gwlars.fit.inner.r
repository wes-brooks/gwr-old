gwlars.fit.inner = function(x, y, coords, indx=NULL, loc, bw=NULL, dist=NULL, s=NULL, mode.select='', tuning=FALSE, predict=FALSE, simulation=FALSE, verbose=FALSE, gwr.weights=NULL, prior.weights=NULL, gweight=NULL, longlat=FALSE, adapt=FALSE, interact, precondition, N=1, tau=3, shrunk.fit, AICc) {
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
        x.interacted = cbind(x, interacted)
        colnames(x.interacted) = c(oldnames, newnames)
    }

    if (mode.select=='CV') { 
        xx = as.matrix(x[-colocated,])
        if (interact) {xx.interacted = as.matrix(x.interacted[-colocated,])}
        yy = as.matrix(y[-colocated])
        w <- prior.weights[-colocated] * gwr.weights[-colocated]  
    } else {
        xx = as.matrix(x)
        if (interact) {xx.interacted = as.matrix(x.interacted)}
        yy = as.matrix(y)
        w <- prior.weights * gwr.weights
    }
    
    if (sum(gwr.weights)==length(colocated)) { return(list(loss.local=Inf, resid=Inf)) } 

    n <- nrow(xx)
    weighted = which(w>0)
    n.weighted = length(weighted)
    
    xx = xx[weighted,]
    if (interact) {xx.interacted = xx.interacted[weighted,]}
    yy = as.matrix(yy[weighted])
    w = w[weighted]

    int.list = list()
    coef.list = list()
    coef.unshrunk.list=list()   
    coef.unshrunk.interacted.list=list()    
    ssr.local = NA

    for (i in 1:N) {
        #Final permutation is the original ordering of the data:
        if (i==N) {            
            permutation = 1:n.weighted
        } else {
            permutation = sample(1:n.weighted, replace=TRUE)
        }

        colocated = which(gwr.weights[weighted][permutation]==1)
        sqrt.w <- diag(sqrt(w[permutation]))        
        yyy = sqrt.w %*% yy[permutation,]
        meany = sum((w*yy)[permutation])/sum(w[permutation])
        yyy = yyy #- meany   
        #normy = sqrt(sum(yyy**2))
        #yyy = yyy / normy
     
        xxx = sqrt.w %*% xx[permutation,]
        if (interact) {xxx.interacted = sqrt.w %*% xx.interacted[permutation,]}

        if (precondition==TRUE) {
            s = svd(xxx)
            F = s$u  %*% diag(1/sqrt(s$d**2 + tau))  %*%  t(s$u)
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
            
            lm.step = try(lm(yyy~xs)) #-1))  # mle fit on standardized
        
            if(class(lm.step) == "try-error") { 
                cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
                return(return(list(loss.local=Inf, resid=Inf)))
            }
            
            beta.lm = lm.step$coeff[-1]                   # mle except for intercept
            adapt.weight = abs(beta.lm)               # weights for adaptive lasso
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(adapt.weight[k])) {
                    xs[,k] = xs[,k] * adapt.weight[k]
                } else {
                    xs[,k] = rep(0, dim(xs)[1])
                    adapt.weight[k] = 0 #This should allow the lambda-finding step to work.
                }
            }
            predx = as.matrix((xx[permutation,] - meanx) * adapt.weight / normx)
                        
            if (interact) {
				one <- rep(1, nrow(xxx.interacted))
				meanx.interacted <- drop(one %*% xxx.interacted) / nrow(xxx.interacted)
				x.interacted.centered <- scale(xxx.interacted, meanx.interacted, FALSE)         # first subtracts mean
				normx.interacted <- sqrt(drop(one %*% (x.interacted.centered**2)))
				names(normx.interacted) <- NULL
				xs.interacted = x.interacted.centered
				
				for (k in 1:dim(x.interacted.centered)[2]) {
					if (normx.interacted[k]!=0) {
						xs.interacted[,k] = xs.interacted[,k] / normx.interacted[k]
					} else {
						xs.interacted[,k] = rep(0, dim(xs.interacted)[1])
						normx.interacted[k] = Inf #This should allow the lambda-finding step to work.
					}
				}
			
				lm.step = try(lm(yyy~xs.interacted)) #-1))  # mle fit on standardized
		
				if(class(lm.step) == "try-error") { 
					cat(paste("Couldn't make a model for finding the SSR at location ", i, ", bandwidth ", bw, "\n", sep=""))
					return(return(list(loss.local=Inf, resid=Inf)))
				}
			
				beta.lm = lm.step$coeff[-1]                   # mle except for intercept
				adapt.weight.interacted = abs(beta.lm)               # weights for adaptive lasso
				for (k in 1:dim(x.interacted.centered)[2]) {
					if (!is.na(adapt.weight.interacted[k])) {
						xs.interacted[,k] = xs.interacted[,k] * adapt.weight.interacted[k]
					} else {
						xs.interacted[,k] = rep(0, dim(xs.interacted)[1])
						adapt.weight.interacted[k] = 0 #This should allow the lambda-finding step to work.
					}
				}
				predx.interacted = as.matrix((xx.interacted[permutation,] - meanx.interacted) * adapt.weight.interacted / normx.interacted)
			}
            
        } else {
            meanx = rep(0, ncol(x))
            adapt.weight = rep(1, ncol(x))
            normx = rep(1, ncol(x))
    
            xs=xxx
            predx = xx[permutation,]
        }
    
        fitx = cbind(sqrt(w[permutation]), xs)
        fity = yyy
    
        model = lars(x=fitx, y=fity, type='lar', normalize=FALSE, intercept=FALSE)
        nsteps = length(model$lambda) + 1   
    
        if (mode.select=='CV') {
        	reps = length(colocated)
        	predx = t(apply(matrix(xx[colocated,], nrow=reps, ncol=dim(xx)[2]), 1, function(X) {(X-meanx) * adapt.weight / normx}))
            vars = apply(predict(model, type='coef')[['coefficients']], 1, function(x) {which(abs(x)>0)})
            df = sapply(vars, length) + 1                        

            predictions = predict(model, newx=predx, type='fit', mode='step')[['fit']]
            loss = colSums(abs(matrix(predictions - matrix(y[colocated], nrow=reps, ncol=nsteps), nrow=reps, ncol=nsteps)))                
            s2 = sum(w[permutation]*lsfit(y=predy, x=predx, wt=w[permutation])$residuals**2) / sum(w[permutation]) 

            loss.local = loss   
            k = which.min(loss)     

        } else if (mode.select=='AIC') {
            #predx = cbind(1, t(apply(xx[permutation,], 1, function(X) {(X-meanx) * adapt.weight / normx})))
            predx = cbind(1, as.matrix(xx[permutation,]))
            predy = as.matrix(yy[permutation])

            if (sum(w[permutation]) > ncol(x)) {               
                #Get the un-penalized intercept
                wcm = apply(predx[,-1], 2, function(x) {sum(x*w[permutation])})/sum(w[permutation])
                coefs = as.matrix(coef(model))
                coefs = t(apply(coefs, 1, function(x) {x * c(1, adapt.weight) * c(1, 1/normx)}))
                coefs[,1] = apply(coefs, 1, function(x) {x[1] - sum(x[-1] * wcm)})
                fitted = predx %*% t(coefs)
                coefs[,1] = coefs[,1] + apply(fitted, 2, function(x) {sum(w[permutation]*(predy-x)) / sum(w[permutation])})
                fitted = predx %*% t(coefs)    
				vars = apply(coefs, 1, function(x) {which(abs(x[-1])>0)})
				df = sapply(vars, length) + 1
                
                s2 = sum(w[permutation]*(fitted[,ncol(fitted)] - yy[permutation])**2) / (sum(w[permutation]) - ncol(x) - 1)#ncol(x))  
                loss = as.vector(apply(fitted, 2, function(z) {sum(w[permutation]*(z - yy[permutation])**2)})/s2 + log(s2) + 2*df)
                k = which.min(loss)
                fitted = fitted[,k]
                localfit = fitted[colocated]      
                df = df[k]                         

                if (length(vars[[k]]) > 0) {
                    varset = vars[[k]] #- 1
                    modeldata = data.frame(y=yy[permutation], xx[permutation,varset])
                    m = lm(y~., data=modeldata, weights=w[permutation])                    
                    result = tryCatch({
                        Xh = diag(sqrt(w[permutation])) %*% as.matrix(cbind(rep(1,length(permutation)), xx[permutation,varset]))
                        H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
                        Hii = H[colocated,colocated]
                    }, error = function(e) {
                        Hii = nrow(x) - 2
                    })
                    
                    if (!shrunk.fit) {
                        fitted = m$fitted
                        localfit = fitted[colocated]
                        df = length(varset) + 1
                        s2 = sum((m$residuals*w[permutation])**2) / (sum(w[permutation]) - df)     
                    }

                    coefs.unshrunk = rep(0, ncol(xx) + 1)
                    coefs.unshrunk[c(1, varset + 1)] = coef(m)
                    s2.unshrunk = sum(m$residuals**2)/(sum(w[permutation]) - length(varset) - 1)

                    se.unshrunk = rep(0, ncol(x) + 1)
                    se.unshrunk[c(1, varset + 1)] = summary(m)$coefficients[,'Std. Error']
                    
                    if (interact) {
						varset.interacted = vars[[k]] #- 1
						for (j in vars[[k]]) {
							varset.interacted = c(varset.interacted, ncol(x)+2*(j-1)+1, ncol(x)+2*j)
						}			

						modeldata = data.frame(y=yy[permutation], xx.interacted[permutation,varset.interacted])
						m = lm(y~., data=modeldata, weights=w[permutation])
						result = tryCatch({
                            Xh = diag(sqrt(w[permutation])) %*% as.matrix(cbind(rep(1,length(permutation)), xx.interacted[permutation,varset.interacted]))
                            H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
                            Hii = H[colocated,colocated]
                        }, error = function(e) {
                            Hii = nrow(x) - 2
                        })                
                        if (!shrunk.fit) {
                            fitted = m$fitted
                            localfit = fitted[colocated]
                            df = length(varset.interacted) + 1
                            s2 = sum((m$residuals*w[permutation])**2) / (sum(w[permutation]) - df)                            
                        }

						coefs.unshrunk.interacted = rep(0, ncol(xx.interacted) + 1)
						coefs.unshrunk.interacted[c(1, varset.interacted + 1)] = coef(m)
						s2.unshrunk.interacted = sum(m$residuals**2)/(sum(w[permutation]) - length(varset.interacted) - 1)

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
                    
                    s2.unshrunk = sum(fity**2)/(sum(w[permutation]) - 1)
                    se.unshrunk = rep(0, ncol(xx) + 1)
                    se.unshrunk[1] = sqrt(s2.unshrunk)
                    
                    coefs.unshrunk.interacted = c(meany, rep(0, ncol(xx.interacted)))
                    se.unshrunk.interacted = c(se.unshrunk[1], rep(0, ncol(xx.interacted)))
                    s2.unshrunk.interacted = s2.unshrunk
                    Hii = 1 / sum(w[permutation])
                }
                
                if (length(colocated)>0) {
                    if (!AICc) {loss.local = sum((w[permutation]*(fitted - yy[permutation])**2)[colocated])/s2 + log(s2) + 2*df/sum(w[permutation])}
                    else {
                        loss.local = Hii
                        ssr.local = sum((w[permutation]*(fitted - yy[permutation])**2)[colocated])
                    }
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

        } else if (mode.select=='BIC') {
            #predx = cbind(1, t(apply(xx[permutation,], 1, function(X) {(X-meanx) * adapt.weight / normx})))
            predx = cbind(1, as.matrix(xx[permutation,]))
            predy = as.matrix(yy[permutation])

            if (sum(w[permutation]) > ncol(x)) {               
                #Get the un-penalized intercept
                wcm = apply(predx[,-1], 2, function(x) {sum(x*w[permutation])})/sum(w[permutation])
                coefs = as.matrix(coef(model))
                coefs = t(apply(coefs, 1, function(x) {x * c(1, adapt.weight) * c(1, 1/normx)}))
                coefs[,1] = apply(coefs, 1, function(x) {x[1] - sum(x[-1] * wcm)})
                fitted = predx %*% t(coefs)
                coefs[,1] = coefs[,1] + apply(fitted, 2, function(x) {sum(w[permutation]*(predy-x)) / sum(w[permutation])})
                fitted = predx %*% t(coefs)    
				vars = apply(coefs, 1, function(x) {which(abs(x[-1])>0)})
				df = sapply(vars, length) + 1
                
                s2 = sum(w[permutation]*(fitted[,ncol(fitted)] - yy[permutation])**2) / (sum(w[permutation]) - ncol(x) - 1)  
                loss = as.vector(apply(fitted, 2, function(z) {sum(w[permutation]*(z - yy[permutation])**2)})/s2 + log(s2) + log(sum(w[permutation]))*df)
                k = which.min(loss)
                fitted = fitted[,k]
                localfit = fitted[colocated]      
                df = df[k]                         

                if (length(vars[[k]]) > 0) {
                    varset = vars[[k]] #- 1
                    modeldata = data.frame(y=yy[permutation], xx[permutation,varset])
                    m = lm(y~., data=modeldata, weights=w[permutation])
                    result = tryCatch({
                        Xh = diag(sqrt(w[permutation])) %*% as.matrix(cbind(rep(1,length(permutation)), xx[permutation,varset]))
                        H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
                        Hii = H[colocated,colocated]
                    }, error = function(e) {
                        Hii = nrow(x) - 2
                    })
                    if (!shrunk.fit) {
                        fitted = m$fitted
                        localfit = fitted[colocated]
                        df = length(varset) + 1
                        s2 = sum((m$residuals*w[permutation])**2) / (sum(w[permutation]) - df)     
                    }

                    coefs.unshrunk = rep(0, ncol(xx) + 1)
                    coefs.unshrunk[c(1, varset + 1)] = coef(m)
                    s2.unshrunk = sum(m$residuals**2)/(sum(w[permutation]) - length(varset) - 1)

                    se.unshrunk = rep(0, ncol(x) + 1)
                    se.unshrunk[c(1, varset + 1)] = summary(m)$coefficients[,'Std. Error']
                    
                    if (interact) {
						varset.interacted = vars[[k]] #- 1
						for (j in vars[[k]]) {
							varset.interacted = c(varset.interacted, ncol(x)+2*(j-1)+1, ncol(x)+2*j)
						}			

						modeldata = data.frame(y=yy[permutation], xx.interacted[permutation,varset.interacted])
						m = lm(y~., data=modeldata, weights=w[permutation])
						
                        result = tryCatch({
                            Xh = diag(sqrt(w[permutation])) %*% as.matrix(cbind(rep(1,length(permutation)), xx.interacted[permutation,varset.interacted]))
                            H = Xh %*% solve(t(Xh) %*% Xh) %*% t(Xh)
                            Hii = H[colocated,colocated]
                        }, error = function(e) {
                            Hii = nrow(x) - 2
                        })              
                        if (!shrunk.fit) {
                            fitted = m$fitted
                            localfit = fitted[colocated]
                            df = length(varset.interacted) + 1
                            s2 = sum((m$residuals*w[permutation])**2) / (sum(w[permutation]) - df)                            
                        }

						coefs.unshrunk.interacted = rep(0, ncol(xx.interacted) + 1)
						coefs.unshrunk.interacted[c(1, varset.interacted + 1)] = coef(m)
						s2.unshrunk.interacted = sum(m$residuals**2)/(sum(w[permutation]) - length(varset.interacted) - 1)

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
                    
                    s2.unshrunk = sum(fity**2)/(sum(w[permutation]) - 1)
                    se.unshrunk = rep(0, ncol(xx) + 1)
                    se.unshrunk[1] = sqrt(s2.unshrunk)
                    
                    coefs.unshrunk.interacted = c(meany, rep(0, ncol(xx.interacted)))
                    se.unshrunk.interacted = c(se.unshrunk[1], rep(0, ncol(xx.interacted)))
                    s2.unshrunk.interacted = s2.unshrunk
                    Hii = 1 / sum(w[permutation])
                }
                
                if (length(colocated)>0) {
                    if (!AICc) {loss.local = sum((w[permutation]*(fitted - yy[permutation])**2)[colocated])/s2 + log(s2) + log(sum(w[permutation]))*df/sum(w[permutation])}
                    else {
                        loss.local = Hii
                        ssr.local = sum((w[permutation]*(fitted - yy[permutation])**2)[colocated])
                    }
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
    }
    
    if (tuning) {
        return(list(loss.local=loss.local, ssr.local=ssr.local, s=s.optimal, sigma2=s2, nonzero=colnames(x)[vars[[s.optimal]]], weightsum=sum(w)))
    } else if (predict) {
        return(list(loss.local=loss.local, coef=coefs))
    } else if (simulation) {
        return(list(loss.local=loss.local, coef=coefs, coeflist=coef.list, s=s.optimal, bw=bw, sigma2=s2, coef.unshrunk=coefs.unshrunk, coef.unshrunk.interacted=coefs.unshrunk.interacted, s2.unshrunk=s2.unshrunk, s2.unshrunk.interacted=s2.unshrunk.interacted, coef.unshrunk.list=coef.unshrunk.list, coef.unshrunk.interacted.list=coef.unshrunk.interacted.list, se.unshrunk=se.unshrunk, se.unshrunk.interacted=se.unshrunk.interacted, nonzero=colnames(x)[vars[[s.optimal]]], fitted=localfit, weightsum=sum(w), loss=loss))
    } else {
        return(list(model=model, loss=loss, coef=coefs, coef.unshrunk=coefs.unshrunk, coeflist=coef.list, s=s.optimal, loc=loc, bw=bw, meanx=meanx, meany=meany, coef.scale=adapt.weight/normx, df=df, loss.local=loss.local, sigma2=s2, sum.weights=sum(w), N=N, fitted=localfit, coef.unshrunk.interacted.list=coef.unshrunk.interacted.list, s2.unshrunk.interacted=s2.unshrunk.interacted, coef.unshrunk.interacted=coefs.unshrunk.interacted, se.unshrunk.interacted=se.unshrunk.interacted))
    }
}
