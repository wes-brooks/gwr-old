#Import the plotting functions:
library(plotrix)
library(sp)
library(fossil)
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


neighbor.weight = function(q, D=NULL, this.coords=NULL, obs.coords=NULL, longlat=FALSE, weight.function, verbose=FALSE, tol=.Machine$double.eps^0.25) {
    if (is.null(D)) {
        bbox <- cbind(range(obs.coords[, 1]), range(obs.coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0

        this.coords = as.numeric(this.coords)
    
        if (!longlat) {
            nrow = dim(obs.coords)[1]
            D2 = (matrix(as.numeric(this.coords), nrow, 2, byrow=TRUE) - obs.coords)**2
            D2 = apply(D2, 1, sum)
            D = sqrt(D2)
        } else {
            n = dim(obs.coords)[1]
            D = sapply(1:dim(obs.coords)[1], function(x) {deg.dist(this.coords[1], this.coords[2], obs.coords[x,1], obs.coords[x,2])})
        }
    }

    beta1 <- min(D)
    beta2 <- 10*max(D)

    optimize(neighbor.diff, lower=beta1, upper=beta2, maximum=FALSE, tol=tol, 
                D=D, weight.function=weight.function, q.target=q, verbose=verbose)
}

neighbor.diff = function(bw, D, q.target, weight.function, verbose) {
    if (verbose) {cat(paste("bandwidth: ", bw, '\n', sep=''))}

    nobs = length(D)
    W = weight.function(D, bw)
    q.this = sum(W)/nobs

    if (verbose) {cat(paste("total weight: ", q.this, '\n', sep=''))}

    return(abs(q.target - q.this))
}


gwr.heatmap <- function(model, variable=NULL, type='coef', cs1=c(0,1), cs2=c(0,1), cs3=c(0,1)) { 
    #Prepare something for plotting:
    if (type=='coef') {
        name.var = variable
        out = vector()
    
        if (is.null(model[['coef.scale']])) {
            coefs=sapply(1:length(model[['s']]), function(k) {coef(sim.model[['model']][[k]])[,which(sim.model[['s.range']]==sim.model[['s']][k])]})
            out = c(out, coefs[name.var,])
        } else {
            coefs=sapply(1:length(model[['s']]), function(k) {model[['coef.scale']][[i]][[name.var]] * coef(sim.model[['model']][[k]])[,which(sim.model[['s.range']]==sim.model[['s']][k])]})
            out = c(out, coefs[name.var,])
        }

    } else if (type=='errors') {
        out = sapply(1:length(model[['s']]), function(k) {t(as.matrix(as.numeric(model[['model']][[k]]$beta[,which(model[['model']][[k]]$lambda==model[['s']][k])]))) %*% as.matrix(as.numeric(c(1, model[['data']][k,c('X1', 'X2', 'Z')]),nrow=4,ncol=1)) + model[['model']][[k]][['a0']][which(model[['model']][[k]]$lambda==model[['s']][k])]}) - eta
    } else if (type=='residuals') {
        out = model[['data']][,model[['response']]] - sapply(1:length(model[['s']]), function(k) {t(as.matrix(as.numeric(model[['model']][[k]]$beta[,which(model[['model']][[k]]$lambda==model[['s']][k])]))) %*% as.matrix(as.numeric(c(1, model[['data']][k,c('X1', 'X2', 'Z')]),nrow=4,ncol=1)) + model[['model']][[k]][['a0']][which(model[['model']][[k]]$lambda==model[['s']][k])]})
    } else if (type=='fit') {
        out = sapply(1:length(model[['s']]), function(k) {t(as.matrix(as.numeric(model[['model']][[k]]$beta[,which(model[['model']][[k]]$lambda==model[['s']][k])]))) %*% as.matrix(as.numeric(c(1, model[['data']][k,c('X1', 'X2', 'Z')]),nrow=4,ncol=1)) + model[['model']][[k]][['a0']][which(model[['model']][[k]]$lambda==model[['s']][k])]})
    } else if (type=='data') {
        out = model[['data']][,model[['response']]]
    }
    
    df.plot = data.frame(output=out)

    #return(df.plot)

    #Isolate the variable to plot:
    locations = model[['coords']]
    coef.surface = as.data.frame(cbind(locations, out))
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

    gwr.matplot(mat, cs1, cs2, cs3, border=NA, show.legend=TRUE, yrev=FALSE, axes=TRUE, ann=TRUE)
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


plot.coef.gwglmnet = function(model, var, breaks=NULL) {
    #Prepare something for plotting:
    name.var = var
    var = vector()
    col.out = which(names(data)==model[['response']])
    
    for (i in 1:length(model[['model']])) {
        m = predict(model[['model']][[i]], s=model[['s']][[i]], newx=data[,-col.out], type='coefficients')
        test = abs(m@x[which(m@i==(which(m@Dimnames[[1]]==name.var)-1))])
        if(name.var %in% m@Dimnames[[1]][m@i+1]) {
            var = c(var, m@x[which(m@i==(which(m@Dimnames[[1]]==name.var)-1))])
        } else {
            var = c(var, 0)
        }
    }
   
    df.plot = data.frame(output=var)
              
    #Put the county names into a form that can be matched.
    locs = unique(merge(model[['data']][,c('x','y','county','state')], model[['coords']]))

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

