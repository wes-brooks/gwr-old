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