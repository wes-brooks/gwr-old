library(ggplot2)
library(sp)
library(scales)
library(grid)

#Bootstrap the source_https function
library(devtools)
install_github("wesesque/brooks")

#Import the data
brooks::source_https('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty-data.r')

#Prepare something for plotting:
year = 1970
df = pov2[pov2$year==year,]
polygons = map_data('county')

#Some counties have their centroid outside their borders or are grouped into a cluster:
#Pepin county:
df[df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',c('x','y')] = c(-92.1048, 44.5823)

#Shawano county:
cluster_id = which(df$COUNTY=="WI_CLUSTER")
df = rbind(df, df[cluster_id,])
df[nrow(df),c('x','y')] = c(-88.707733, 44.788658)

#Oconto county:
df = rbind(df, df[cluster_id,])
df[nrow(df),c('x','y')] = c(-88.014221, 44.877282)

locs = data.frame(df[,c('x','y')])
locs = unique(locs)

#Plotting constants:
col.bg='gray85'
response = "logitindpov"
group = 'group'
borderlines=NULL
col.outline='white'


#Merge the polygons with the locs:
mergedata = data.frame()
for (k in unique(polygons[,group])) {
    id = which(point.in.polygon(pol.x=polygons[polygons[,group]==k,1], pol.y=polygons[polygons[,group]==k,2], point.x=locs[,1], point.y=locs[,2])==1)
    if (length(id)==1) {
        shape = polygons[polygons[,group]==k,]
        shape = cbind(shape, id, df[id,])
        mergedata = rbind(mergedata, shape)
    }
}


#Draw the map of poverty rate
ll = range(mergedata[,response], na.rm=TRUE)
map <- ggplot(mergedata, aes(long,lat,group=group)) +
    geom_polygon(aes_string(fill=response)) +
    scale_fill_gradient(high='orange', low='white', limits=ll, na.value='gray50', name="") +
    coord_map(project='globular') +
    theme(panel.background=element_rect(fill=col.bg, colour=col.outline)) +
    ggtitle(column.map[[response]]) +
    theme(plot.margin=unit(c(0,0,0,0), "cm"), legend.margin=unit(0,'cm')) +
    scale_x_continuous('') +
    scale_y_continuous('')



#######################################
#Choropleth maps of the raw covariates:
maps = list()
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof')

for (v in predictors) {
    ll = range(mergedata[,v], na.rm=TRUE)

    #Draw the map of proportion in agriculture
    maps[[v]] <- ggplot(mergedata, aes(long,lat,group=group)) +
        geom_polygon(aes_string(fill=v)) +
        scale_fill_gradient(high='orange', low='white', limits=ll, na.value='gray50', name="") +
        coord_map(project='globular') +
        theme(panel.background=element_rect(fill=col.bg, colour=col.outline)) +
        ggtitle(column.map[[v]]) +
        theme(plot.margin=unit(c(0,0,0,0), "cm"), legend.margin=unit(0, "cm")) +
        scale_x_continuous('') +
        scale_y_continuous('')
}

pdf("~/git/gwr/figures/practice-talk/poverty-response.pdf", width=6, height=6)
map
dev.off()

pdf("~/git/gwr/figures/practice-talk/poverty-covariates.pdf", width=11, height=6)
multiplot(plotlist=maps, cols=3)
dev.off()