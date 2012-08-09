library(spgwr)
library(lars)
library(glmnet)
library(maps)
library(ggplot2)
library(fossil)
library(lqa)


#extract reference data
mapcounties <- map_data("county")
mapstates <- map_data("state")


#Import the plotting functions:
setwd("~/git/gwr/code")
source("utils.r")

#Import poverty data
#pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
pov = read.csv("~/git/gwr/data/dropbox/data.csv", header=TRUE)
pov$y1960 = 1 - pov$y1970 - pov$y1980 - pov$y1990 - pov$y2000 - pov$y2006
years = c('y1960', 'y1970', 'y1980', 'y1990', 'y2000', 'y2006')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty', pov_denom='population', poor='number of poor',
    logitfampov='logit( proportion families in poverty)', pfire='pfire')



#Process the poverty data so that each column appears only once and the year is added as a column.
pov2 = list()
for (column.name in names(column.map)) {
    pov2[[column.name]] = pov[[column.name]]
}

#Find the columns we haven't yet matched:
"%w/o%" <- function(x, y) x[!x %in% y]
missed = names(pov) %w/o% outer(names(column.map), years, FUN=function(x, y) {paste(x, y, sep="")})

for (column.name in missed) {
    col = pov[,column.name]
    pov2[[column.name]] = col
}

#Add the year column to the pov2 data list.
pov2[['year']] = vector()
for (i in 1:dim(pov)[1]) {
    if (pov$y1960[i]) {
        y = 1960
    } else if (pov$y1970[i]) {
        y=1970
    } else if (pov$y1980[i]) {
        y=1980
    } else if (pov$y1990[i]) {
        y=1990
    } else if (pov$y2000[i]) {
        y=2000
    } else if (pov$y2006[i]) {
        y=2006
    } 
    pov2[['year']] = c(pov2[['year']], y)
}

#Convert pov2 from a list to a data frame:
pov2 = data.frame(pov2)

#Set the matrix that we'll use for GWR
#df = pov2[pov2$year==2006,]
df = pov2


#Define which variables we'll use as predictors of poverty:
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro', 'pind')
f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))
f2 = as.formula(paste("pindpov ~ ", paste(predictors, collapse="+"), sep=""))

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}

model.data = pov2[,predictors]
model.data[['logitindpov']] = pov2[['logitindpov']]
#model.data[['population']] = pov2[['pov_denom']]

#Put the locations into a more friendly format:
loc = data.frame(x=df$x, y=df$y, county=as.character(df$county), state=as.character(df$state))
loc.unique = unique(loc)
model.locs = loc.unique[,c('county', 'state')]
model.coords = loc.unique[,c('x', 'y')]

#Select the bandwidth for a GWR model without selection:
#k.nn = gwr.sel(f, data=df, coords=loc[,c('x'0,'y')], longlat=TRUE, adapt=TRUE, gweight=gwr.bisquare)
k.nn = 0.09446181
#bw = gwr.sel(f, data=df, coords=loc[,c('x','y')], longlat=TRUE, adapt=FALSE, gweight=gwr.bisquare)
#bw = 295.9431  #bandwidth in kilometers
#bw = 125.314 #from gwlars.sel
bw = 127 #-ish, from gw.adalars


#Use the spgwr package to produce a model without selection
knn_model = gwr(f, data=df, coords=as.matrix(df[,c('x','y')]), adapt=k.nn, longlat=TRUE, gweight=gwr.bisquare, fit.points=as.matrix(model.coords))
bw_model = gwr(f, data=df, coords=as.matrix(df[,c('x','y')]), bandwidth=bw, longlat=TRUE, gweight=gwr.bisquare, fit.points=as.matrix(model.coords))

#Plot the full model on a map
plot.coef.gwr(model=bw_model, var='pind', locs=loc)









#Compute the matrix of distances (in kilometers)
n = dim(loc)[1]
D = as.matrix(earth.dist(loc[,c('x','y')]),n,n)
w = bisquare(D, bw=bw)

#Do the same for just the unique locations
n.unique = dim(model.coords)[1]
D.unique = as.matrix(earth.dist(model.coords),n,n)
w.unique = bisquare(D.unique, bw=bw)


cv_error = data.frame()
w.lasso.geo = list()
glmnet.geo = list()
adaglmnet.geo = list()
lqa.geo = list()
adalars.geo = list()
coefs = list()
ss = seq(0, 1, length.out=100)
lambda = seq(0, 5, length.out=5000)
l.lars = vector()
l.glmnet = vector()
l.adaglmnet = vector()
l.lqa = vector()
l.adalars = vector()
col.out = which(names(model.data)=='logitindpov')
reps = dim(model.data)[1]/n.unique


for(i in 1:dim(model.coords)[1]) {
    colocated = which(loc$x==model.coords$x[i] & loc$y==model.coords$y[i])
    loow = w.unique[i,-colocated]

    #model = lm(f, data=model.data[-colocated,], weights=loow)
    
    #Weight is based on the combination of the distances and the county populations.
    w.sqrt <- diag(rep(sqrt(loow), reps) * sqrt(df$pov_denom[-colocated]))
    #block.eig = list()
    #block.eig[['vectors']] = matrix(0, length(loow)*reps, length(loow)*reps)
    #block.eig[['vectors']][,1:length(loow)] = rbind(w.eig$vectors, w.eig$vectors, w.eig$vectors, w.eig$vectors, w.eig$vectors, w.eig$vectors)
    #block.eig[['values']] = c(reps*w.eig$values, rep(0, (reps-1)*length(loow)))

    #w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% t(w.eig$vectors)
    #w.sqrt = block.eig[['vectors']] %*% diag(sqrt(block.eig[['values']])) %*% t(block.eig[['vectors']])
    
    w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[-colocated,predictors]), y=w.sqrt %*% as.matrix(df$logitindpov[-colocated]))
    glmnet.geo[[i]] = glmnet(x=as.matrix(df[-colocated, predictors]), y=as.matrix(cbind(df$pindpov[-colocated], 1-df$pindpov[-colocated])), weights=rep(loow, reps) * df$pov_denom[-colocated], family='binomial')
    #lqa.geo[[i]] = lqa(x=as.matrix(df[-colocated, predictors]), y=as.matrix(cbind(df$pindpov[-colocated], 1-df$pindpov[-colocated])), weights=rep(loow, reps), family='binomial', penalty=adaptive.lasso)
    
    #Compute the adaptive lasso estimates
    #Generate the adaptive lasso weights
    x.weighted = w.sqrt %*% as.matrix(model.data[-colocated,predictors])
    y.weighted = w.sqrt %*% as.matrix(model.data$logitindpov[-colocated])
    m<-ncol(x.weighted)
    n<-nrow(x.weighted)
    one <- rep(1, n)
    meanx <- drop(one %*% x.weighted)/n
    x.centered <- scale(x.weighted, meanx, FALSE)         # first subtracts mean
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

    #xs <- scale(x.centered, FALSE, normx)        # now rescales with norm (not sd)
    
    out.ls = lm(y.weighted~xs)                      # ols fit on standardized
    beta.ols = out.ls$coeff[2:(m+1)]       # ols except for intercept
    ada.weight = abs(beta.ols)                      # weights for adaptive lasso

    for (k in 1:dim(x.centered)[2]) {
        if (!is.na(ada.weight[k])) {
            xs[,k] = xs[,k] * ada.weight[k]
        } else {
            xs[,k] = rep(0, dim(xs)[1])
            ada.weight[k] = 0 #This should allow the lambda-finding step to work.
        }
    }

    for (k in 1:dim(x.centered)[2]) {
        if (!is.na(ada.weight[k])) {
            xs[,k] = xs[,k] * ada.weight[k]
        } else {
            xs[,k] = rep(0, dim(xs)[1])
            ada.weight[k] = 0 #This should allow the lambda-finding step to work.
        }
    }
    
    #Use the lars algorithm to fit the model
    adalars.geo[[i]] = lars(x=xs, y=y.weighted, normalize=FALSE)


    ############################
    #Adaptive lasso for glmnet:
    #############################
    #Compute the adaptive lasso estimates
    #Generate the adaptive lasso weights
    m<-ncol(x)
    n<-nrow(x)
    one <- rep(1, n)
    meanx <- drop(one %*% x)/n
    x.centered <- scale(x, meanx, FALSE)         # first subtracts mean
    normx <- sqrt(drop(one %*% (x.centered^2)))
    names(normx) <- NULL
    xs.adaglmnet = x.centered
    for (k in 1:dim(x.centered)[2]) {
        if (normx[k]!=0) {
            xs.adaglmnet[,k] = xs.adaglmnet[,k] / normx[k]
        } else {
            xs.adaglmnet[,k] = rep(0, dim(xs.adaglmnet)[1])
            normx[k] = Inf #This should allow the lambda-finding step to work.
        }
    }

    df.adaglmnet = data.frame(y=y,xs.adaglmnet)
    
    if (!is.null(weights)) {
        rawglm = glm(y~., data=df.adaglmnet[-colocated,], family='binomial', weights=loow*weights[-colocated])
    } else {
        rawglm = glm(y~., data=df.adaglmnet[-colocated,], family='binomial')
    }
    beta.glm = rawglm$coeff[2:(m+1)]                      # ols except for intercept
    adaglmnet.weight = abs(beta.glm)                      # weights for adaptive lasso
    # weights for adaptive lasso
    for (k in 1:dim(x.centered)[2]) {
        if (!is.na(adaglmnet.weight[k])) {
            xs.adaglmnet[,k] = xs.adaglmnet[,k] * adaglmnet.weight[k]
        } else {
            xs.adaglmnet[,k] = rep(0, dim(xs.adaglmnet)[1])
            adaglmnet.weight[k] = 0 #This should allow the lambda-finding step to work.
        }
    }

    df = data.frame(y=y,xs.adaglmnet)


    adaglmnet.geo[[i]] = glmnet(x=as.matrix(df[-colocated, predictors]), y=as.matrix(cbind(df$pindpov[-colocated], 1-df$pindpov[-colocated])), weights=rep(loow, reps) * df$pov_denom[-colocated], family='binomial')

    l.lars = c(l.lars, which.min(colSums(abs(predict(w.lasso.geo[[i]], newx=model.data[colocated,-col.out], s=lambda, type='fit', mode='lambda')[['fit']] - model.data[colocated,col.out])))/1000)
    l.glmnet = c(l.glmnet, glmnet.geo[[i]][['lambda']][which.min(colSums(abs(predict(glmnet.geo[[i]], newx=as.matrix(model.data[colocated,-col.out]), type='response') - model.data[colocated,col.out])))])
    l.adaglmnet = c(l.adaglmnet, adaglmnet.geo[[i]][['lambda']][which.min(colSums(abs(predict(glmnet.geo[[i]], newx=as.matrix(model.data[colocated,-col.out]), type='response') - model.data[colocated,col.out])))])
    l.adalars = c(l.adalars, which.min(colSums(abs(predict(adalars.geo[[i]], newx=scale(model.data[colocated,-col.out], center=meanx, scale=normx/ada.weight), type='fit', mode='step')[['fit']] - model.data[colocated,col.out]))))

    #l.lqa = c(l.lqa, 
    print(i)
}




