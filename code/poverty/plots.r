## You should run poverty.r to generate the plots for this script.

#Import the county polygons
library(ggplot2)
library(gridExtra)
library(mapproj)

upper_midwest = c('illinois', 'indiana', 'iowa', 'michigan', 'minnesota', 'wisconsin')
county = map_data('county')
county = county[county$region %in% upper_midwest,]

#Establish lists to hold the plots:
plots = list()
plots[['enet']] = list()
plots[['lasso']] = list()

for (yr in years) {
    ####Plot a choropleth of the results:    
    #Correct the locations of some small counties (only affect plotting)
	cluster_id = which(df$COUNTY=="WI_CLUSTER")
	n.counties = nrow(df)
    year = as.character(yr)

    #for (select in c("lasso", "enet")) {
    for (select in c("enet")) {
        #Pepin county:
        model[[select]][[year]][['coords']][df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',] = c(-92.1048, 44.5823)

        #Shawano county:
        model[[select]][[year]][['model']][['models']][[n.counties+1]] = model[[select]][[year]][['model']][['models']][[cluster_id]]
        model[[select]][[year]][['coords']] = rbind(model[[select]][[year]][['coords']], c(-88.707733, 44.788658))

        #Oconto county:
        model[[select]][[year]][['model']][['models']][[n.counties+2]] = model[[select]][[year]][['model']][['models']][[cluster_id]]
        model[[select]][[year]][['coords']] = rbind(model[[select]][[year]][['coords']], c(-88.014221, 44.877282)) 
    
        plots[[select]][[year]] = list()
        for (v in predictors) {
            plots[[select]][[year]][[v]] = plot.gwselect(model[[select]][[year]], part='coef.unshrunk', var=v, polygons=county, title=v, col.bg='gray85') + theme(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
        }

        #pp = plots[[select]][[year]]
        dev.new()
        #pdf(paste('~/git/gwr/figures/poverty/', yr, '-', select, '-linear-coefficients-unshrunk.pdf', sep=''), width=8, height=16)
        grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], ncol=2)
        #dev.off()
    }
}