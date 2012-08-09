library(spgwr)
data(columbus)

#Import the plotting functions:
setwd("~/git/gwr/code")
source("matplot.r")
source("legend.r")

#Import poverty data
pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty',
    logitfampov='logit( proportion families in poverty)')

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

#Convert pov2 from a list to a data frame.
pov2 = data.frame(pov2)

#Correct the Y2K bug
pov2 = within(pov2, year <- as.numeric(as.character(year)) + 1900)
pov2 = within(pov2, year <- ifelse(year<1960, year+100, year))

#Define the bisquare weight function:
bisquare = function(x, z, bw) {
    ifelse( r(x=x, z=c(z[1], z[2])) < bw, (1 - (r(x=x, z=c(z[1],z[2]))/ bw)**2)**2, 0)
}

bisquare = function(R, bw) {
    ifelse( R < bw, (1 - (R/bw)**2)**2, 0)
}

#Define the distance function
r = function(x, z) { 
    sqrt((x[1]-z[1])**2 + (x[2]-z[2])**2)
}

#Fix an x.i and get the distance to all other x's:
R = function(x.i, xy.mat) {
    R = sapply(1:dim(xy.mat)[1], FUN=function(i) {r(x.i, xy.mat[i,])} )
    return(R)
}

W = function(x.i, xy.mat, bw) {
    distance = R(x.i, xy.mat)
    W = bisquare(distance, bw)
    return(W)
}

#Define a grid of locations where we'll fit a GWR model:
n=20
xx = as.vector(quantile(df$x, 1:n/(n+1)))
yy = as.vector(quantile(df$y, 1:n/(n+1)))
locs = cbind(x=rep(xx,each=n), y=rep(yy,times=n))

#Use this trick to compute the matrix of distances very quickly
n = dim(pov)[1]
D1 = matrix(rep(pov$x,n), n,n)
D2 = matrix(rep(pov$y,n), n,n)
D = sqrt((D1-t(D1))**2 + (D2-t(D2))**2)

#Define which variables we'll use as predictors of poverty:
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}

#Use the lasso for GWR models of poverty with 2006 data:
df = pov2[pov2$year==2006,]
w.lasso.geo = list()
coefs = list()
ss = seq(0, 1, length.out=100)

for(i in 1:dim(df)[1]) {
    w = bisquare(D[,i], bw=3)

    model = lm(f, data=df, weights=w)
    
    w.eig <- eigen(diag(w))
    w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[,predictors]), y=as.matrix(df$logitindpov))
    
    for (col in predictors) {
        coefs[[col]] = c(coefs[[col]], model$coef[[col]])
    }
    
    print(i)
}


#Use the methods of spgwr to select a bandwidth and fit a GWR model for poverty:
bw = gwr.sel(income~1, data=columbus, coords=cbind(x,y), adapt=FALSE, gweight=gwr.bisquare)
gwr.model1 = gwr(income~1, data=columbus, coords=cbind(x,y), bandwidth=bw, gweight=gwr.bisquare, hatmatrix=TRUE)

#Make a GAM model for poverty:
gm.pag = gam(logitindpov ~ s(pag, k=10, by=as.matrix(ll)), data=df)


#Set the directory to store plots
plot_dir = "figures/"

#plot the model outputs:
for(pred in predictors) {    
    #Isolate the variable to plot:
    coef.surface = as.data.frame(cbind(locs[,c("x", "y")], get(pred)))
    names(coef.surface)[3] = pred
    
    #Heatmap of the data
    locations = with(coef.surface, list(lat=unique(y), long=unique(x)))
    mat = matrix(NA, nrow=length(locations[['lat']]), ncol=length(locations[['long']]))
    rownames(mat) <- sort(unique(coef.surface$y), decreasing=F)
    colnames(mat) <- sort(unique(coef.surface$x), decreasing=F)
    
    #Put the coefficients into a lat-long matrix
    for(row in 1:dim(coef.surface)[1]) {
        mat[as.character(coef.surface[row,"y"]), as.character(coef.surface[row,"x"])] = ifelse(!is.na(coef.surface[row,pred]), coef.surface[row,pred], NA)
    }

    #par(bty='n')
    matplot(mat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=TRUE, yrev=FALSE, axes=TRUE, ann=TRUE)
}

