library(spgwr)
library(lars)

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




#Limit the data to just 100 points for now for computational reasons
pov3 = pov2[pov2$year==2006,][1:100,]
df = pov3




#Define a grid of locations where we'll fit a GWR model:
n=20
xx = as.vector(quantile(df$x, 1:n/(n+1)))
yy = as.vector(quantile(df$y, 1:n/(n+1)))
locs = cbind(x=rep(xx,each=n), y=rep(yy,times=n))

#Use this trick to compute the matrix of distances very quickly
n = dim(pov3)[1]
D1 = matrix(rep(pov3$x,n), n,n)
D2 = matrix(rep(pov3$y,n), n,n)
D = sqrt((D1-t(D1))**2 + (D2-t(D2))**2)
w = bisquare(D, bw=3)


#Define which variables we'll use as predictors of poverty:
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}

cv_error = data.frame()
w.lasso.geo = list()
coefs = list()


w.lasso.geo = list()
coefs = list()
ss = seq(0, 1, length.out=100)
lambda = seq(0, 2, length.out=2000)

for(i in 1:dim(df)[1]) {
    w = D[i,-i]

    model = lm(f, data=df[-i,], weights=w)
    
    w.eig <- eigen(diag(w))
    w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[-i,predictors]), y=as.matrix(df$logitindpov[-i]))
    
    for (col in predictors) {
        coefs[[col]] = c(coefs[[col]], model$coef[[col]])
    }
    
    cv_error = rbind(cv_error, df$logitindpov[i] - predict(w.lasso.geo[[i]], newx=df[i,predictors], type='fit', mode='lambda', s=lambda)$fit)

    print(i)
}

names(cv_error) = lambda

er = vector()
for (k in 1:dim(cv_error)[2]) {
    er = c(er, sum((cv_error**2)[k]))
    cat(paste(names(cv_error)[k], ": ", sum((cv_error**2)[k]), "\n", sep=""))
}



for (j in 1:dim(df)[1]) {
    loodf = df[-j,]
    looD = D[-j,-j]

    ss = seq(0, 1, length.out=100)
    w.lasso.geo[[j]] = list()
    coefs[[j]] = list()
    
    affected = which(w[j,] > 0)

    for(i in affected) {    
        model = lm(f, data=loodf, weights=w)
        
        w.eig <- eigen(diag(w))
        w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
        w.lasso.geo[[j]][[i]] = lars(x=w.sqrt %*% as.matrix(loodf[,predictors]), y=as.matrix(loodf$logitindpov))
        
        for (col in predictors) {
            coefs[[j]][[col]] = c(coefs[[j]][[col]], model$coef[[col]])
        }        
        print(i)
    }    
    print(j)
}


#Plot the CV error along the lasso path for each data point
pred = data.frame()
l = seq(0, 1, length.out=1000)

for(k in 1:100) {
    pred = rbind(pred, abs(predict(w.lasso.geo[[k]], newx=pov3[k,predictors], mode='lambda', s=l)$fit - pov3$logitindpov[k]))
}