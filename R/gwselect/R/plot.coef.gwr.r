plot.coef.gwr = function(model, var, locs, breaks=NULL) {
    #Prepare something for plotting:
    locs = merge(unique(locs), model$SDF@coords)
    name.var = var
    var = model$SDF@data[,name.var]
   
    df.plot = data.frame(output=var)
    
    #Put the county names into a form that can be matched.
    df.plot$county = tolower(as.character(locs$county))
    df.plot$state = tolower(as.character(locs$state))
    for (i in 1:dim(df.plot)[1]) {
        county = gsub("['-. ]", '', df.plot$county[i])
        df.plot$county[i] = paste(county, tolower(df.plot$state[i]), sep=',')
    }    
    
    #extract reference data
    mapcounties <- map_data('county')
    mapstates <- map_data('state')
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=','))
    mergedata <- merge(midwestcounties, df.plot, by.x='county', by.y='county')
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    
    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste("Coefficient of '", name.var, "' in a model for logitindpov without LARS", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}
    