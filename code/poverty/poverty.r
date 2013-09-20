#Import external libraries
library(gwselect)
library(spgwr)
registerCores(n=3)
library(gridExtra)
library(mapproj)

#Import the data
source('code/poverty/poverty-data.r')

#Establish lists to hold the bandwidths
bw = list()
bw[['enet']] = list()
bw[['lasso']] = list()

#Establish lists to hold the models
model = list()
model[['enet']] = list()
model[['lasso']] = list()

years = c(1960, 1970, 1980, 1990, 2000, 2006)

for (yr in years) {
    #Use the lasso for GWR models of poverty with 2006 data:
    df = pov2[pov2$year==yr,]
    df = merge(x=df, y=pops, by.x="fips", by.y="FIPS")
    year = as.character(yr)

    #Define which variables we'll use as predictors of poverty:
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
    f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))

    bw[['lasso']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, mode.select="BIC", gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)
    model[['lasso']][[year]] = gwglmnet(formula=f, data=df, family='gaussian', alpha=1, coords=df[,c('x','y')], longlat=TRUE, N=1, mode.select='BIC', bw=bw[['lasso']][[year]], gweight=bisquare, method='dist', simulation=TRUE, adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

    bw[['enet']][[year]] = gwglmnet.sel(formula=f, data=df, family='gaussian', alpha='adaptive', coords=df[,c('x','y')], longlat=TRUE, mode.select="BIC", gweight=bisquare, tol=0.01, s=NULL, method='dist', adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)
    model[['enet']][[year]] = gwglmnet(formula=f, data=df, family='gaussian', alpha='adaptive', coords=df[,c('x','y')], longlat=TRUE, N=1, mode.select='BIC', bw=bw[['enet']][[year]], gweight=bisquare, method='dist', simulation=TRUE, adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

    #f.spgwr = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))
    #bw.spgwr = gwr.sel(formula=f.spgwr, data=df, coords=as.matrix(df[,c('x','y')]), gweight=gwr.bisquare, method="aic", show.error.messages=TRUE)
    #model.spgwr = gwr(formula=f.spgwr, data=df, coords=as.matrix(df[,c('x','y')]), bandwidth=bw.spgwr, gweight=gwr.bisquare)
}

#Make a map using polygons provided to the plot function.              
#Put the county names into a form that can be matched.
county = map_data('county')

plots = list()
plots[['enet']] = list()
plots[['lasso']] = list()

for (yr in years) {
	#Correct the locations of some small counties (only affect plotting)
	cluster_id = which(df$COUNTY=="WI_CLUSTER")
	n.counties = nrow(df)
    year = as.character(yr)

    for (select in c("lasso", "enet")) {
        model[[select]][[year]][['model']][['models']][[n.counties+1]] = model[[select]][[year]][['model']][['models']][[cluster_id]]
        model[[select]][[year]][['coords']][df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',] = c(-92.1048, 44.5823)
        model[[select]][[year]][['coords']] = rbind(model[[select]][[year]][['coords']], c(-88.7285, 44.9291))
    
        plots[[select]][[year]] = list()
        for (v in predictors) {
            plots[[select]][[year]][[v]] = plot.gwselect(model[[select]][[year]], part='coef.unshrunk', var=v, polygons=county, title=v, col.bg='gray85') + theme(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
        }

        pp = plots[[select]][[year]]
        dev.new()
        pdf(paste('figures/poverty/', year, '.', select, '.linear.coefficients.unshrunk.pdf', sep=''), width=8, height=16)
        grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], pp[['pwh']], pp[['pblk']], pp[['phisp']], pp[['metro']], ncol=2)
        dev.off()
    }
}