## You should run poverty.r to generate the plots for this script.

#Import the county polygons
library(ggplot2)
library(gridExtra)
library(mapproj)
county = map_data('county')

#Establish lists to hold the plots:
plots = list()
plots[['enet']] = list()
plots[['lasso']] = list()

for (yr in years) {
    ####Plot a choropleth of the results:    
    #Correct the locations of some small counties (only affect plotting)
	cluster_id = which(df$COUNTY=="WI_CLUSTER")
	n.counties = nrow(df)

    #for (select in c("lasso", "enet")) {
    for (select in c("enet")) {
        model[[select]][[yr]][['model']][['models']][[n.counties+1]] = model[[select]][[yr]][['model']][['models']][[cluster_id]]
        model[[select]][[yr]][['coords']][df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',] = c(-92.1048, 44.5823)
        model[[select]][[yr]][['coords']] = rbind(model[[select]][[yr]][['coords']], c(-88.7285, 44.9291))
    
        plots[[select]][[year]] = list()
        for (v in predictors) {
            plots[[select]][[yr]][[v]] = plot.gwselect(model[[select]][[yr]], part='coef.unshrunk', var=v, polygons=county, title=v, col.bg='gray85') + theme(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
        }

        pp = plots[[select]][[year]]
        #dev.new()
        pdf(paste('figures/poverty/', yr, '-', select, '-linear-coefficients-unshrunk.pdf', sep=''), width=8, height=16)
        #grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], pp[['pwh']], pp[['pblk']], pp[['phisp']], pp[['metro']], ncol=2)
        grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], ncol=2)
        dev.off()
    }
}