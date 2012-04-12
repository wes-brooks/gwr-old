#Import the plotting functions:
library(plotrix)
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
