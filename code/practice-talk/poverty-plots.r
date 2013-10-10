library(ggplot2)
library(sp)
library(scales)
library(grid)

#Import the data
setwd('~/git/gwr/code/poverty')
source('poverty-data.r')
source("~/git/brooks/code/multiplot.r")

#Prepare something for plotting:
df = pov2[pov2$year==1970,]
polygons = map_data('county')

#Some counties have their centroid outside their borders or are grouped into a cluster
cluster_id = which(df$COUNTY=="WI_CLUSTER")
n.counties = nrow(df)
df = rbind(df, df[cluster_id,])
df[df$STATE=='Wisconsin' & df$COUNTY=='PEPIN',c('x','y')] = c(-92.1048, 44.5823)
df[nrow(df),c('x','y')] = c(-88.7285, 44.9291)

locs = data.frame(df[,c('x','y')])
locs = unique(locs)

#Plotting constants:
col.bg='gray85'
name.var = "logitindpov"
group = 'group'
title = 'Logit of county poverty rate, 1970'
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
name.var = 'logitindpov'
map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes_string(fill=name.var))
map <- map + scale_fill_gradient(high='orange', low='white', limits=range(mergedata[,name.var], na.rm=TRUE), na.value='gray50', name="") + coord_map(project='globular')
map <- map + theme(panel.background=element_rect(fill=col.bg, colour=col.outline))

#Annotate the map with borderlines
if (!is.null(borderlines)) {map <- map + geom_path(data=borderlines, colour='white', size=0.75)}

#Plot the map
map <- map + ggtitle(title) + theme(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')


name.var = 'pag'
#Draw the map of proportion in agriculture
map2 <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes_string(fill=name.var))
map2 <- map2 + scale_fill_gradient(high='orange', low='white', limits=range(mergedata[,name.var], na.rm=TRUE), na.value='gray50', name="") + coord_map(project='globular')
map2 <- map2 + theme(panel.background=element_rect(fill=col.bg, colour=col.outline))

#Annotate the map with borderlines
if (!is.null(borderlines)) {map2 <- map2 + geom_path(data=borderlines, colour='white', size=0.75)}

#Plot the map
map2 <- map2 + ggtitle("Proportion working in agriculture") + theme(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')



name.var = 'pman'
#Draw the map of proportion in agriculture
map3 <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes_string(fill=name.var))
map3 <- map3 + scale_fill_gradient(high='orange', low='white', limits=range(mergedata[,name.var], na.rm=TRUE), na.value='gray50', name="") + coord_map(project='globular')
map3 <- map3 + theme(panel.background=element_rect(fill=col.bg, colour=col.outline))

#Annotate the map with borderlines
if (!is.null(borderlines)) {map2 <- map2 + geom_path(data=borderlines, colour='white', size=0.75)}

#Plot the map
map3 <- map3 + ggtitle("Proportion working in manufacturing") + theme(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')


jpeg("../../figures/practice-talk/poverty-data.jpg", width=12, height=4, units='in', res=72)
multiplot(map, map2, map3, cols=3)
dev.off()