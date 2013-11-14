#Import external libraries
require(MASS)
require(gwselect)
require(spgwr)
require(devtools)
registerCores(n=3)

#Import the data
source_url('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty-data.r')

yr = 1970
year = as.character(yr)
df = pov2[pov2$year==yr,]

####Produce the models via lasso and via elastic net:
#Define which variables we'll use as predictors of poverty:
#predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')
f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))



resampled = gwglmnet.ridge(formula=f, data=df, family='gaussian', bw=bw[['GWAL']][['1970']], coords=df[,c('x','y')], longlat=TRUE, gweight=bisquare, tol=0.01, parallel=TRUE, interact=TRUE, verbose=TRUE)
