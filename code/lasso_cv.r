library(spgwr)
library(lars)
library(maps)
library(ggplot2)
library(fossil)


#extract reference data
mapcounties <- map_data("county")
mapstates <- map_data("state")


#Import the plotting functions:
setwd("~/git/gwr/code")
source("utils.r")

#Import poverty data
pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty',
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

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}

model.data = pov2[,predictors]
model.data[['logitindpov']] = pov2[['logitindpov']]

#Put the locations into a more friendly format:
loc = data.frame(x=df$x, y=df$y)
loc.unique = unique(loc)
counties = data.frame(x=df$x, y=df$y, county=as.character(df$COUNTY), state=as.character(df$STATE))

#Select the bandwidth for a GWR model without selection:
#k.nn = gwr.sel(f, data=df, coords=loc, longlat=TRUE, adapt=TRUE, gweight=gwr.bisquare)
k.nn = 0.09446181
#bw = gwr.sel(f, data=df, coords=loc, longlat=TRUE, adapt=FALSE, gweight=gwr.bisquare)
bw = 295.9431  #bandwidth in kilometers

#Use the spgwr package to produce a model without selection
knn_model = gwr(f, data=df, coords=loc, adapt=k.nn, longlat=TRUE, gweight=gwr.bisquare, fit.points=unique(loc))
bw_model = gwr(f, data=df, coords=loc, bandwidth=bw, longlat=TRUE, gweight=gwr.bisquare, fit.points=unique(loc))

#Plot the full model on a map
plot.coef.gwr(model=bw_model, var='pind')









#Compute the matrix of distances (in kilometers)
n = dim(loc)[1]
D = as.matrix(earth.dist(loc),n,n)
w = bisquare(D, bw=bw)

#Do the same for just the unique locations
n.unique = dim(loc.unique)[1]
D.unique = as.matrix(earth.dist(loc.unique),n,n)
w.unique = bisquare(D.unique, bw=bw)


cv_error = data.frame()
w.lasso.geo = list()
coefs = list()
ss = seq(0, 1, length.out=100)
lambda = seq(0, 2, length.out=2000)
l = vector()
col.out = which(names(model.data)=='logitindpov')
reps = dim(loc)[1]/dim(loc.unique)[1]


for(i in 1:dim(loc.unique)[1]) {
    colocated = which(loc$x==loc.unique$x[i] & loc$y==loc.unique$y[i])
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
     
    l = c(l, which.min(colSums(abs(predict(w.lasso.geo[[i]], newx=model.data[colocated,-col.out], s=lambda, type='fit', mode='lambda')[['fit']] - model.data[colocated,col.out])))/1000)
    print(i)
}


