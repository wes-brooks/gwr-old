library(gwselect)
registerCores(n=7)

#Import poverty data
pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
pov$X09pop = as.numeric(gsub(",", "", as.character(pov$X09pop)))
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty',
    logitfampov='logit( proportion families in poverty)', pfire='proportion financial, insurance, real estate')

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

#Use the lasso for GWR models of poverty with 2006 data:
df = pov2[pov2$year==2006,]

#Define which variables we'll use as predictors of poverty:
#weights=rep(1, nrow(pov2))
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("logitindpov ~ -1+ ", paste(predictors, collapse="+"), sep=""))
#bw.pov = gwlars.sel(formula=f, data=df, coords=df[,c('x','y')], weights=df$X09pop, longlat=TRUE, gweight=bisquare, mode='step', mode.select='AIC', s=NULL, method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)
bw.pov.noweights = gwlars.sel(formula=f, data=df, coords=df[,c('x','y')], longlat=TRUE, gweight=bisquare, mode='step', mode.select='BIC', s=NULL, method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)
#bw.pov (pov2, CV, knn) = 0.189408989047768
#bw.pov (pov2, AIC, knn) = 0.491
#bw.pov (df, AIC, knn) Bandwidth: 0.866496634540107. Loss: 0.116838909368888
model.pov = gwlars(formula=f, data=df, coords=df[,c('x','y')], longlat=TRUE, gweight=bisquare, bw=bw.pov.noweights, mode='step', mode.select='BIC', s=NULL, method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)



#Define which variables we'll use as predictors of poverty:
#weights=rep(1, nrow(pov2))
#predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
#f = as.formula(paste("pindpov ~ ", paste(predictors, collapse="+"), sep=""))
#bw = gwglmnet.sel(formula=f, data=pov2, family='binomial', weights=weights, coords=pov2[,c('x','y')], adapt=FALSE, gweight=bisquare, s=NULL, method="knn", tol=0.001, longlat=TRUE, parallel=FALSE, verbose=FALSE, precondition=FALSE)
#bw = 0.007443209

#Define which variables we'll use as predictors of poverty:
weights=rep(1, nrow(pov2))
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("pindpov ~ -1 +", paste(predictors, collapse="+"), sep=""))
bw = 0.007443209
model = gwglmnet(formula=f, data=pov2, family='binomial', weights=weights, bw=bw, coords=pov2[,c('x','y')], adapt=TRUE, gweight=bisquare, s=NULL, method="knn", tol=0.001, longlat=TRUE, parallel=FALSE, verbose=FALSE, precondition=FALSE)



#Make a map using polygons provided to the plot function.
          
#Put the county names into a form that can be matched.
county = map_data('county')
#state = map_data('state')

plot.gwselect(model.pov, var="pfire", polygons=county)
