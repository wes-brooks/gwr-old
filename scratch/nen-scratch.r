#Import external libraries
require(devtools)
install_github('gwselect', 'wesesque')
require(gwselect)
require(RCurl)
registerCores(n=3)

#Import the data
source_url('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty-data.r')

#Establish lists to hold the bandwidths
bw = list()
bw[['GWEN']] = list()
bw[['GWAL']] = list()
bw[['GWR']] = list()

#Establish lists to hold the models
model = list()
model[['GWEN']] = list()
model[['GWAL']] = list()
model[['GWR']] = list()

years = c(1960, 1970, 1980, 1990, 2000, 2006)
years = c(1970)


yr=1970


    #Isolate one year of data
    year = as.character(yr)
    df = pov2[pov2$year==yr,]
    
    ####Produce the models via lasso and via elastic net:
    #Define which variables we'll use as predictors of poverty:
    #predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')
    f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))



bw[['GWAL']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, mode.select="BIC", gweight=spherical, tol=1, s=NULL, method='dist', adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE, AICc=TRUE)
