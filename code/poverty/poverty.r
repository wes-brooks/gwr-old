#Import external libraries
require(gwselect)
require(spgwr)
require(devtools)
registerCores(n=3)

#If the 'brooks' package isnt loaded then import it from github:
if (!'package:brooks' %in% search()) { install_github('wesesque/brooks') }

#Import the data
brooks::source_https('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty-data.r')

#Establish lists to hold the bandwidths
bw = list()
bw[['GWEN']] = list()
bw[['GWAL']] = list()

#Establish lists to hold the models
model = list()
model[['GWEN']] = list()
model[['GWAL']] = list()

years = c(1960, 1970, 1980, 1990, 2000, 2006)
years = c(1970)

for (yr in years) {
    #Isolate one year of data
    year = as.character(yr)
    df = pov2[pov2$year==yr,]
    
    ####Produce the models via lasso and via elastic net:
    #Define which variables we'll use as predictors of poverty:
    #predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')
    f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))

    #Lasso model
    #bw[['GWAL']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, mode.select="BIC", gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)
    #model[['GWAL']][[year]] = gwglmnet(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, N=1, mode.select='BIC', bw=bw[['lasso']][[year]], gweight=bisquare, method='dist', simulation=TRUE, adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

    #Elastic net model:
    bw[['GWEN']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', alpha='adaptive', coords=df[,c('x','y')], longlat=TRUE, mode.select="BIC", gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)
    model[['GWEN']][[year]] = gwglmnet(formula=f, data=df, family='gaussian', alpha='adaptive', coords=df[,c('x','y')], longlat=TRUE, N=1, mode.select='BIC', bw=bw[['enet']][[year]], gweight=bisquare, method='dist', simulation=TRUE, adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

    #f.spgwr = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))
    #bw.spgwr = gwr.sel(formula=f.spgwr, data=df, longlat=TRUE, coords=as.matrix(df[,c('x','y')]), gweight=gwr.bisquare, method="aic", show.error.messages=TRUE)
    #model.spgwr = gwr(formula=f.spgwr, data=df, longlat=TRUE, coords=as.matrix(df[,c('x','y')]), bandwidth=bw.spgwr, gweight=gwr.bisquare)

    #Use my code to do the traditional GWR; (currently buggy)
    #oracle = lapply(1:533, function(x) {return(predictors)})
    #bw.gwr = gwglmnet.sel(f, data=df, oracle=oracle, coords=df[,c('x','y')], family='gaussian', alpha=1, mode.select='BIC', longlat=TRUE, gweight=bisquare, tol=0.01, method='dist', parallel=FALSE, interact=FALSE, verbose=TRUE, shrunk.fit=FALSE, AICc=TRUE)
    #model.gwr = gwglmnet(f, data=df, oracle=oracle, coords=df[,c('x','y')], family='gaussian', alpha=1, mode.select='BIC', longlat=TRUE, bw=bw.gwr, gweight=bisquare, method='dist', simulation=TRUE, parallel=FALSE, interact=FALSE, verbose=TRUE, shrunk.fit=FALSE, AICc=TRUE)
}