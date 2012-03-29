library(spgwr)
data(columbus)

pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
heads = c('pindpov', 'logitindpov', 'pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'pind', 'phisp', 'metro', 'pfampov', 'logitfampov')
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty',
    logitfampov='logit( proportion families in poverty)')

pov2 = list()

for (column.name in heads) {
    col = vector()

    for (year in years) {
        if (paste(column.name, year, sep="") %in% names(pov)) {
            indx = which(names(pov)==paste(column.name, year, sep=""))
            col = c(col, pov[,indx])
        }
        else {
            col = c(col, rep(NA, dim(pov)[1]))
        }
    }
    pov2[[column.name]] = col
}

#Find the columns we haven't yet matched:
"%w/o%" <- function(x, y) x[!x %in% y]
missed = names(pov) %w/o% outer(heads, years, FUN=function(x, y) {paste(x, y, sep="")})

for (column.name in missed) {
    col = rep(pov[,column.name], length(years))
    pov2[[column.name]] = col
}

pov2[['year']] = vector()
for (year in years) {
    pov2[['year']] = c(pov2[['year']], rep(year, dim(pov)[1]))
}

pov2 = data.frame(pov2)

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


predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))

for (col in predictors) {
    assign(col, vector())
}

df = pov2[pov2$year==2006,]

for(i in 1:dim(df)[1]) {
    w = W(df[i, c('x', 'y')], df[,c('x','y')], bw=1)

    model = lm(f, data=df, weights=w)
    
    w.eig <- eigen(diag(w))
    w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    w.lasso[[i]] = lars(x=w.sqrt %*% as.matrix(df[,predictors]), y=as.matrix(df$logitindpov))

    for (col in predictors) {
        coef = get(col)
        assign(col, c(coef, model$coef[[col]]))
    }
    
    print(i)
}

bw = gwr.sel(income~1, data=columbus, coords=cbind(x,y), adapt=FALSE, gweight=gwr.bisquare)
gwr.model1 = gwr(income~1, data=columbus, coords=cbind(x,y), bandwidth=bw, gweight=gwr.bisquare, hatmatrix=TRUE)

