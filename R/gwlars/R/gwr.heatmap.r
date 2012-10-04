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