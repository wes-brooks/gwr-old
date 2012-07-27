library(maptools)
wi=readShapePoly("~/git/gwr/data/northernWisconsin/plss_quarter-sections/qscppoly")
ashqs = which(wi@data$DTRSQ %in% landcover$ID)
ash = wi@polygons[ashqs]

ash.coords = data.frame()
for (k in 1:length(ashqs)) {
    ash.coords = rbind(ash.coords, cbind(91.24667 - atan((ash[[k]]@Polygons[[1]]@coords[,1]-423538)/4402443), 46.35056 + atan((ash[[k]]@Polygons[[1]]@coords[,2]-653773)/6378100), rep(wi@data$DTRSQ[ashqs[k]],nrow(ash[[k]]@Polygons[[1]]@coords))))
}

#aa = ash.coords
#aa$x = -aa$x
#names(aa) = c('x','y','DTRSQ')

#mergedata <- merge(aa, landcover, by.x='DTRSQ', by.y='ID')

#map <- ggplot(mergedata, aes(x,y,group=DTRSQ)) + geom_polygon(aes(fill=out))
#map <- map + scale_fill_gradient(low='white', high='red', limits=range(mergedata$out, na.rm=TRUE), name='coef') + coord_map(project='globular')