#Import the plotting functions:
library(plotrix)
library(sp)
setwd('~/git/gwr/code')
source('matplot.r')
#source('legend.r')

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

#Fix an x.i and get the distance to all other x’s:
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
    #Prepare something for plotting:
    name.var = variable
    var = vector()
    
    for (i in 1:length(model[['model']])) {
        var = c(var, model[['coef.scale']][i][[1]][[name.var]] * coef.lars(model[['model']][[i]], mode=model[['mode']], s=model[['s']][i])[[name.var]])
    }
   
    df.plot = data.frame(output=var)


    #Isolate the variable to plot:
    locations = model[['coords']]
    coef.surface = as.data.frame(cbind(locations, var))
    names(coef.surface)[3] = name.var
    
    #Heatmap of the data
    locations = list(lat=unique(locations[,2]), long=unique(locations[,1]))
    mat = matrix(NA, nrow=length(locations[['lat']]), ncol=length(locations[['long']]))
    rownames(mat) <- sort(unique(coef.surface[,2]), decreasing=F)
    colnames(mat) <- sort(unique(coef.surface[,1]), decreasing=F)         
    
    #Put the coefficients into a lat-long matrix
    for(row in 1:dim(coef.surface)[1]) {
        mat[as.character(coef.surface[row,2]), as.character(coef.surface[row,1])] = 
            ifelse(!is.na(coef.surface[row,name.var]), coef.surface[row,name.var], NA)
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
        mat[as.character(coef.surface[row,'y']), as.character(coef.surface[row,'x'])] = 
            ifelse(!is.na(coef.surface[row,variable]), coef.surface[row,variable], NA)
    }

    #par(bty='n')
    gwr.matplot(mat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=TRUE, yrev=FALSE, axes=TRUE, ann=TRUE)
}




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
    




plot.coef.gwlars = function(model, var, locs, data, s=NULL, breaks=NULL) {
    #Prepare something for plotting:
    name.var = var
    var = vector()
    
    for (i in 1:length(model)) {
        var = c(var, coef.lars(model[['model']][[i]], mode=model[['mode']], s=model[['s']][i])[[name.var]])
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
    mapcounties <- map_data('county')
    mapstates <- map_data('state')
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    if (mean(df.plot$output, na.rm=TRUE)<=0) {       
        map <- map + scale_fill_gradient(low='red', high='white', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    } else {
        map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    }

    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste("Coefficient of '", name.var, "' in a model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}
    

plot.effect.gwlars = function(model, var, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    name.var = var
    var = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        var = c(var, data[[name.var]][i]*coef.lars(model[[i]], newx=data[i,-col.out], mode='lambda', s=l[i])[[name.var]])
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
    mapcounties <- map_data('county')
    mapstates <- map_data('state')
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=','))
    mergedata <- merge(midwestcounties, df.plot, by.x = 'county', by.y = 'county')
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name='effect') + coord_map(project="globular")
    
    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste(c('Coefficient of ', name.var, ' in a model for logitindpov'), collapse='')) #+ guides(fill=guide_legend(reverse=TRUE))
}
    

plot.residuals.gwlars = function(model, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    err = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        err = c(err, predict.lars(model[[i]], newx=data[i,-col.out], mode='lambda', s=l[i])[['fit']] - data[i,col.out])
    }
   
    df.plot = data.frame(output=err)
              
    #Put the county names into a form that can be matched.
    df.plot$county = tolower(as.character(locs$county))
    df.plot$state = tolower(as.character(locs$state))
    for (i in 1:dim(df.plot)[1]) {
        county = gsub("['-. ]", '', df.plot$county[i])
        df.plot$county[i] = paste(c(county, tolower(df.plot$state[i])), collapse=',')
    }    
    
    #extract reference data
    mapcounties <- map_data('county')
    mapstates <- map_data('state')
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    if (mean(df.plot$output, na.rm=TRUE)<=0)       
        map <- map + scale_fill_gradient(low='red', high='white', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    
    else
        map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    

    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste("Residuals in a GW-LARS model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}


plot.sign.resid.gwlars = function(model, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    err = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        err = c(err, predict.lars(model[[i]], newx=data[i,-col.out], mode='lambda', s=l[i])[['fit']] - data[i,col.out])
    }
   
    df.plot = data.frame(output=sign(err))
              
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
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    #map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    map <- map + scale_fill_gradient(low='red', high='white', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste("Residuals in a GW-LARS model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}



plot.residuals.gwr = function(model, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    err = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        err = c(err, predict.gwr(model[[i]], newx=data[i,-col.out], mode='lambda', s=l[i])[['fit']] - data[i,col.out])
    }
   
    df.plot = data.frame(output=err)
              
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
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    if (mean(df.plot$output, na.rm=TRUE)<=0)      
        map <- map + scale_fill_gradient(low='red', high='white', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    
    else
        map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    

    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste("Residuals in a GW-LARS model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}


plot.coef.gwglmnet = function(model, var, locs, l, data, breaks=NULL) {
    #Prepare something for plotting:
    name.var = var
    var = vector()
    col.out = which(names(data)=='logitindpov')
    
    for (i in 1:length(model)) {
        m = predict(model[[i]], s=l[[i]], newx=data[,-col.out], type='coefficients')
        test = abs(m@x[which(m@i==(which(m@Dimnames[[1]]==name.var)-1))])
        if(name.var %in% m@Dimnames[[1]][m@i+1]) {
            var = c(var, abs(m@x[which(m@i==(which(m@Dimnames[[1]]==name.var)-1))]))
        } else {
            var = c(var, 0)
        }
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
    mapcounties <- map_data('county')
    mapstates <- map_data('state')
    
    #limit our view to the midwest:
    midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.plot$state),]
    midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.plot$state),]
    
    #merge data with ggplot county coordinates
    midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep=","))
    mergedata <- merge(midwestcounties, df.plot, by.x = "county", by.y = "county")
    mergedata <- mergedata[order(mergedata$group, mergedata$order),]
    
    #draw map
    map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
    if (mean(df.plot$output, na.rm=TRUE)<=0)        
        map <- map + scale_fill_gradient(low='red', high='white', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    else
        map <- map + scale_fill_gradient(low='white', high='red', limits=range(df.plot$output, na.rm=TRUE), name='coef') + coord_map(project='globular')
    

    map <- map + opts(panel.background=theme_rect(fill='green', colour='red'))
    
    #add state borders
    map <- map + geom_path(data=midweststates, colour='white', size=0.75)
    
    #add county borders
    map <- map + geom_path(data=midwestcounties, colour='white', size=0.5, alpha=0.1)
    map + opts(title=paste("Coefficient of '", name.var, "' in a model for logitindpov", sep='')) #+ guides(fill=guide_legend(reverse=TRUE))
}

