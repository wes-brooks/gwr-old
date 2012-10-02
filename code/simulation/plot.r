heatdata = matrix(NA, 30,30)
vars = c('X1', 'X2', 'X3', 'X4', 'X5', 'Z')
#pdf("../../figures/simulation/adapt-coefs.pdf", width=6, height=9)
layout(matrix(1:6, 3, 2))

for (var in vars) {
    for (j in 1:length(model[['model']][['models']])) {
        m = model[['model']][['models']][[j]]    
        heatdata[j] = m[['coef']][var,]
        #heatdata = rbind(heatdata, c(m[['loc']], out=m[['coef']][var,]))
    }
    gwr.matplot(heatdata, c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
    title(paste("Coef of ", var, sep=""))
}
#dev.off()

coefdata = matrix(NA, 40,40)
var = 'B1'

for (j in 1:length(model[['model']][['models']])) {
    m = model[['model']][['models']][[j]]    
    heatdata[j] = m[['coef']][var,]
    #heatdata = rbind(heatdata, c(m[['loc']], out=m[['coef']][var,]))
}

matplot(B1, c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
matplot(B2, c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F)



coef.selected = matrix(NA, 30,30)
vars = c('X1', 'X2', 'X3', 'X4', 'X5', 'Z')
#pdf("../../figures/simulation/adapt-coefs.pdf", width=6, height=9)
layout(matrix(1:6, 3, 2))

for (var in vars) {
    for (j in 1:length(model[['model']][['models']])) {
        m = model[['model']][['models']][[j]]    
        coef.selected[j] = ifelse(abs(m[['coef']][var,])>0, 1, 0)
        #heatdata = rbind(coef.selected, c(m[['loc']], out=m[['coef']][var,]))
    }
    gwr.matplot(coef.selected, c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
    title(paste("Coef of ", var, sep=""))
}
#dev.off()