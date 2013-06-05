library(gwselect)
registerCores(n=3)
library(gridExtra)
library(mapproj)

#Import and manipulate data for modeling
source("poverty-data.r")

#Get just the 1970 data:
year = 1970
df = pov2[pov2$year==year,]
df = merge(x=df, y=pops, by.x="fips", by.y="FIPS")   

#Define which variables we'll use as predictors of poverty:
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))

#Make the gwlars model
bw = gwlars.sel(formula=f, data=df, coords=df[,c('x','y')], longlat=TRUE, gweight=bisquare, mode.select='AIC', method="knn", tol=0.001, precondition=FALSE, adapt=TRUE, verbose=FALSE, parallel=TRUE, interact=TRUE, shrunk.fit=FALSE)
model = gwlars(formula=f, data=df, N=200, coords=df[,c('x','y')], longlat=TRUE, gweight=bisquare, bw=bw, mode.select='AIC', s=NULL, method="knn", tol=0.001, precondition=FALSE, adapt=TRUE, verbose=FALSE, parallel=TRUE, interact=TRUE, shrunk.fit=FALSE)

locs=as.matrix(df[,c('x','y')])
coefs = t(sapply(1:dim(df)[1], function(y) {as.vector(model[['model']][['models']][[y]][['coef']])}))
colnames(coefs) = c("(Intercept)", predictors)

#Deal with some counties that are oddly named or have their center of mass outside their borders
cluster_id = which(df$COUNTY=="WI_CLUSTER")
n.counties = nrow(df)
coefs = rbind(coefs, coefs[cluster_id,], coefs[cluster_id,]) #Add entries for Shawano and Oconto counties
locs[df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',] = c(-92.1048, 44.5823)
locs = rbind(locs, c(-88.62121, 44.77598)) #Shawano, WI
locs = rbind(locs, c(-87.8741, 44.887012)) #Oconto, WI
write.csv(locs, "~/git/gwr/output/poverty/centroids.csv", row.names=FALSE)

#Extract the bootstrap estimates of the coefficients, and write them to csv files.
write.csv(coefs, paste("~/git/gwr/output/poverty/coefs-", year, ".csv", sep=""), row.names=FALSE)
bootstrap = list()
for (p in colnames(coefs)) {
    indx = which(colnames(coefs)==p)
    bootstrap[[p]] = t(sapply(1:dim(df)[1], function(y) {sapply(model[['model']][['models']][[y]][['coef.unshrunk.interacted.list']], function(x) {x[indx]})}))
    bootstrap[[p]] = rbind(bootstrap[[p]], bootstrap[[p]][cluster_id,])
    write.csv(bootstrap[[p]], paste("~/git/gwr/output/poverty/coef-bootstrap-", p, "-", year, ".csv", sep=""), row.names=FALSE)
}

#Make a map using polygons provided to the plot function.              
#Put the county names into a form that can be matched.
#county = map_data('county')

#plots[[as.character(year)]] = list()
#plots.unshrunk[[as.character(year)]] = list()
#for (v in predictors) {
#    plots[[as.character(year)]][[v]] = plot.gwselect(model.spgwr, locs=locs, part='coef', values=coefs[[v]], polygons=county, title=v, col.bg='gray85') + opts(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
#}


#pp = plots[[as.character(year)]]
#dev.new()
#pdf(paste('figures/poverty/', year, '.gwr.pdf', sep=''), width=8, height=16)
#grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], pp[['pwh']], pp[['pblk']], pp[['phisp']], pp[['metro']], ncol=2)
#dev.off()


#v = 'pex'
#plot(x=coefs[,v], y=sapply(model[['1970']][['model']][['models']], function(x) {x[['coef']][v,]}), xlab='gwr', ylab='enet')
#abline(a=0,b=1)
