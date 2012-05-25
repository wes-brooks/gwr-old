#Import the plotting functions:
library(plotrix)
library(sp)
setwd("~/git/gwr/code")
source("matplot.r")
#source("legend.r")

#Define the bisquare weight function:
bisquare = function(x, z, bw) {
    ifelse( r(x=x, z=c(z[1], z[2])) < bw, (1 - (r(x=x, z=c(z[1],z[2]))/ bw)**2)**2, 0)
}

bisquare = function(R, bw) {
    ifelse( R < bw, (1 - (R/bw)**2)**2, 0)
}

#Define the distance function
r = function(x, z) { 
    sqrt((x[1]-z[1])**2 + (x[2]-z[2])**2)
}

#Fix an x.i and get the distance to all other x's:
R = function(x.i, xy.mat) {
    R = sapply(1:dim(xy.mat)[1], FUN=function(i) {r(x.i, xy.mat[i,])} )
    return(R)
}

W = function(x.i, xy.mat, bw) {
    distance = R(x.i, xy.mat)
    W = bisquare(distance, bw)
    return(W)
}

gwr.heatmap <- function(model, variable) { 
    #Isolate the variable to plot:
    locations = model$SDF@coords
    coef.surface = as.data.frame(cbind(locations, model$SDF@data[[variable]]))
    names(coef.surface)[3] = variable
    
    #Heatmap of the data
    locations = with(coef.surface, list(lat=unique(y), long=unique(x)))
    mat = matrix(NA, nrow=length(locations[['lat']]), ncol=length(locations[['long']]))
    rownames(mat) <- sort(unique(coef.surface$y), decreasing=F)
    colnames(mat) <- sort(unique(coef.surface$x), decreasing=F)         
    
    #Put the coefficients into a lat-long matrix
    for(row in 1:dim(coef.surface)[1]) {
        mat[as.character(coef.surface[row,"y"]), as.character(coef.surface[row,"x"])] = 
            ifelse(!is.na(coef.surface[row,variable]), coef.surface[row,variable], NA)
    }

    #par(bty='n')
    gwr.matplot(mat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=TRUE, yrev=FALSE, axes=TRUE, ann=TRUE)
}

gwr.heatmap.homebrew <- function(model, variable) {   
    #Isolate the variable to plot:
    locations = model[['coords']]
    coef.surface = as.data.frame(cbind(locations, model[['coefs']][[variable]]))
    names(coef.surface)[3] = variable
    
    #Heatmap of the data
    locations = with(coef.surface, list(lat=unique(y), long=unique(x)))
    mat = matrix(NA, nrow=length(locations[['lat']]), ncol=length(locations[['long']]))
    rownames(mat) <- sort(unique(coef.surface$y), decreasing=F)
    colnames(mat) <- sort(unique(coef.surface$x), decreasing=F)         
    
    #Put the coefficients into a lat-long matrix
    for(row in 1:dim(coef.surface)[1]) {
        mat[as.character(coef.surface[row,"y"]), as.character(coef.surface[row,"x"])] = 
            ifelse(!is.na(coef.surface[row,variable]), coef.surface[row,variable], NA)
    }

    #par(bty='n')
    gwr.matplot(mat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=TRUE, yrev=FALSE, axes=TRUE, ann=TRUE)
}




plot.coef.gwr = function(model, var, locs, breaks=NULL) {
    #Prepare something for plotting:
    locs = merge(unique(locs), full_model$SDF@coords)
    name.var = var
    var = model$SDF@data[,var]
   
    df.plot = data.frame(output=var)
              
    #Put each cell into a color bucket    
    if (is.null(breaks)) {
        breaks = as.vector(quantile(output, (1:6)/6))
    }
    
    bin = vector()
    for (i in 1:dim(df.plot)[1]) {
        bin = c(bin, min(which(breaks >= df.plot$output[i])))
    }

    colors = c("#980043", "#C994C7", "#D4B9DA", "#DD1C77", "#DF65B0", "#F1EEF6")
    df.plot$color = colors[bin]
     

    #Put the county names into a form that can be matched.
    df.plot$county = tolower(as.character(locs$county))
    df.plot$state = tolower(as.character(locs$state))
    for (i in 1:dim(df.plot)[1]) {
        county = gsub("['-. ]", '', df.plot$county[i])
        df.plot$county[i] = paste(county, tolower(df.plot$state[i]), sep=',')
    }
    
    
    #extract reference data
    mapcounties <- map_data("county")
    mapstates <- map_data("state")
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name="coef") + #scale_fill_brewer(palette="PuRd") +
        coord_map(project="globular")
    
    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour="white", size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour="white", size=0.5, alpha=0.1)
    map + opts(title=paste("Coefficient of '", name.var, "' in a model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}
    




plot.coef.gwlars = function(model, var, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    name.var = var
    var = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        var = c(var, coef.lars(model[[i]], newx=data[i,-col.out], mode='lambda', s=1000*l[i])[[name.var]])
    }
   
    df.plot = data.frame(output=var)
              
    #Put the county names into a form that can be matched.
    df.plot$county = tolower(as.character(locs$county))
    df.plot$state = tolower(as.character(locs$state))
    for (i in 1:dim(df.plot)[1]) {
        county = gsub("['-. ]", '', df.plot$county[i])
        df.plot$county[i] = paste(county, tolower(df.plot$state[i]), sep=',')
    }    
    
    #extract reference data
    mapcounties <- map_data("county")
    mapstates <- map_data("state")
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name="coef") + #scale_fill_brewer(palette="PuRd") +
        coord_map(project="globular")
    
    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour="white", size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour="white", size=0.5, alpha=0.1)
    map + opts(title=paste("Coefficient of '", name.var, "' in a model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}
    

plot.effect.gwlars = function(model, var, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    name.var = var
    var = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        var = c(var, data[[name.var]][i]*coef.lars(model[[i]], newx=data[i,-col.out], mode='lambda', s=1000*l[i])[[name.var]])
    }

    df.plot = data.frame(output=var)

    #Put the county names into a form that can be matched.
    df.plot$county = tolower(as.character(locs$county))
    df.plot$state = tolower(as.character(locs$state))
    for (i in 1:dim(df.plot)[1]) {
        county = gsub("['-. ]", '', df.plot$county[i])
        df.plot$county[i] = paste(county, tolower(df.plot$state[i]), sep=',')
    }
    
    #extract reference data
    mapcounties <- map_data("county")
    mapstates <- map_data("state")
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name="coef") + #scale_fill_brewer(palette="PuRd") +
        coord_map(project="globular")
    
    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour="white", size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour="white", size=0.5, alpha=0.1)
    map + opts(title=paste("Coefficient of '", name.var, "' in a model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}
    
