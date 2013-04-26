library(spgwr)
library(gwselect)
library(gridExtra)
library(mapproj)

f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))
bw.spgwr = gwr.sel(formula=f, data=df, coords=as.matrix(df[,c('x','y')]), gweight=gwr.Gauss, method="aic", longlat=TRUE)
model.spgwr = gwr(formula=f, data=df, coords=as.matrix(df[,c('x','y')]), bandwidth=bw.spgwr, gweight=gwr.Gauss, longlat=TRUE)


locs=as.matrix(df[,c('x','y')])
coefs = model.spgwr$SDF@data
#plot.gwselect(model=model.spgwr, values=coefs[['pex']], part='coef', locs=as.matrix(df[,c('x','y')]), polygons=county, col.bg='gray80')


cluster_id = which(df$COUNTY=="WI_CLUSTER")
n.counties = nrow(df)
coefs = rbind(coefs, coefs[cluster_id,])
locs[df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',] = c(-92.1048, 44.5823)
locs = rbind(locs, c(-88.7285, 44.9291))

#Make a map using polygons provided to the plot function.              
#Put the county names into a form that can be matched.
county = map_data('county')

plots[[as.character(year)]] = list()
plots.unshrunk[[as.character(year)]] = list()
for (v in predictors) {
    plots[[as.character(year)]][[v]] = plot.gwselect(model.spgwr, locs=locs, part='coef', values=coefs[[v]], polygons=county, title=v, col.bg='gray85') + opts(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
}


pp = plots[[as.character(year)]]
#dev.new()
pdf(paste('figures/poverty/', year, '.gwr.pdf', sep=''), width=8, height=16)
grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], pp[['pwh']], pp[['pblk']], pp[['phisp']], pp[['metro']], ncol=2)
dev.off()


v = 'pex'
plot(x=coefs[,v], y=sapply(model[['1970']][['model']][['models']], function(x) {x[['coef']][v,]}), xlab='gwr', ylab='enet')
abline(a=0,b=1)
