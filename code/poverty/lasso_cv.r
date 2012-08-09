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
pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
#pov = read.csv("~/git/gwr/data/dropbox/data.csv", header=TRUE)
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty', pov_denom='population', poor='number of poor',
    logitfampov='logit( proportion families in poverty)', pfire='pfire')

#Process the poverty data so that each column appears only once and the year is added as a column.
pov2 = list()
for (column.name in names(column.map)) {
    col = vector()
    for (year in years) {
        if (paste(column.name, year, sep="") %in% names(pov)) {
            indx = which(names(pov)==paste(column.name, year, sep=""))
            col = c(col, pov[,indx])
        }
        else { col = c(col, rep(NA, dim(pov)[1])) }
    }
    pov2[[column.name]] = col
}

#Find the columns we haven't yet matched:
"%w/o%" <- function(x, y) x[!x %in% y]
missed = names(pov) %w/o% outer(names(column.map), years, FUN=function(x, y) {paste(x, y, sep="")})

for (column.name in missed) {
    col = rep(pov[,column.name], length(years))
    pov2[[column.name]] = col
}

#Add the year column to the pov2 data list.
pov2[['year']] = vector()
for (year in years) {
    pov2[['year']] = c(pov2[['year']], rep(year, dim(pov)[1]))
}

#Convert pov2 from a list to a data frame:
pov2 = data.frame(pov2)

#Correct the Y2K bug
pov2 = within(pov2, year <- as.numeric(as.character(year)) + 1900)
pov2 = within(pov2, year <- ifelse(year<1960, year+100, year))



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

#Put the locations into a more friendly format:
loc = data.frame(x=df$x, y=df$y, county=as.character(df$county), state=as.character(df$county))
loc.unique = unique(loc)
model.locs = loc.unique[,c('county', 'state')]
model.coords = loc.unique[,c('x', 'y')]

#Select the bandwidth for a GWR model without selection:
#k.nn = gwr.sel(f, data=df, coords=loc[,c('x','y')], longlat=TRUE, adapt=TRUE, gweight=gwr.bisquare)
k.nn = 0.09446181
#bw = gwr.sel(f, data=df, coords=loc[,c('x','y')], longlat=TRUE, adapt=FALSE, gweight=gwr.bisquare)
#bw = 295.9431  #bandwidth in kilometers
#bw = 125.314 #from gwlars.sel
bw = 127 #-ish, from gw.adalars
#bw.glmnet.adapt = gwglmnet.sel(y~., data=df, coords=loc[,c('x','y')], longlat=TRUE, adapt=TRUE, gweight=gwr.bisquare, weights=df$pov_denom)
bw.glmnet.adapt = 1635.349

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
lqa.geo = list()
adalars.geo = list()
coefs = list()
ss = seq(0, 1, length.out=100)
lambda = seq(0, 5, length.out=5000)
l.lars = vector()
l.glmnet = vector()
l.lqa = vector()
l.adalars = vector()
col.out = which(names(model.data)=='logitindpov')
reps = dim(model.data)[1]/n.unique

adalars.scale = list()
adalars.normx = list()


for(i in 1:dim(model.coords)[1]) {
    colocated = which(loc$x==model.coords$x[i] & loc$y==model.coords$y[i])
    loow = w.unique[i,-colocated]

    #model = lm(f, data=model.data[-colocated,], weights=loow)
    
    w.sqrt <- diag(rep(sqrt(loow), reps))
    #block.eig = list()
    #block.eig[['vectors']] = matrix(0, length(loow)*reps, length(loow)*reps)
    #block.eig[['vectors']][,1:length(loow)] = rbind(w.eig$vectors, w.eig$vectors, w.eig$vectors, w.eig$vectors, w.eig$vectors, w.eig$vectors)
    #block.eig[['values']] = c(reps*w.eig$values, rep(0, (reps-1)*length(loow)))

    #w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% t(w.eig$vectors)
    #w.sqrt = block.eig[['vectors']] %*% diag(sqrt(block.eig[['values']])) %*% t(block.eig[['vectors']])
    
    w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[-colocated,predictors]), y=w.sqrt %*% as.matrix(df$logitindpov[-colocated]))
    glmnet.geo[[i]] = glmnet(x=as.matrix(df[-colocated, predictors]), y=as.matrix(cbind(df$pindpov[-colocated], 1-df$pindpov[-colocated])), weights=rep(loow, reps), family='binomial')
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
    adalars.normx[[i]] = normx
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
    adalars.scale[[i]] = ada.weight
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

    l.lars = c(l.lars, which.min(colSums(abs(matrix(predict(w.lasso.geo[[i]], newx=as.matrix(model.data[colocated,-col.out]), s=lambda, type='fit', mode='lambda')[['fit']] - matrix(model.data[colocated,col.out]), nrow=reps, ncol=length(lambda)))))/1000)
    l.glmnet = c(l.glmnet, glmnet.geo[[i]][['lambda']][which.min(colSums(abs(predict(glmnet.geo[[i]], newx=as.matrix(model.data[colocated,-col.out]), type='response') - model.data[colocated,col.out])))])
    l.adalars = c(l.adalars, lambda[which.min(colSums(abs(matrix(predict(adalars.geo[[i]], newx=scale(as.matrix(model.data[colocated,-col.out], nrow=reps, ncol=dim(model.data[,-col.out])[2]), center=meanx, scale=normx/ada.weight), type='fit', mode='lambda', s=lambda)[['fit']] - model.data[colocated,col.out], nrow=reps, ncol=length(lambda)))))])

    #l.lqa = c(l.lqa, 
    print(i)
}




