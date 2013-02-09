selection.aggregate = list()
vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
settings = 1:18

for (setting in settings) {
    selection = list()

    nsims = ifelse(setting==18,99,100)
    for (k in 1:nsims) {
        sim = (setting-1)*100 + k

        #Import our coefficient estimates
        filename = paste("output/output/CoefEstimates.", cluster, ".", sim, ".csv", sep="")
        estimates = read.csv(filename, header=TRUE)
        colnames(estimates) = vars

        for (v in vars) {
            #Calculate how often each variable was selected for inclusion in the local models
            col = which(colnames(estimates) == v)
            if (k==1) {
                selection[[v]] = as.matrix(ifelse(estimates[,col]==0, 0, 1))
            } else {
                selection[[v]] = cbind(selection[[v]], as.matrix(ifelse(estimates[,col]==0, 0, 1)))
            }
        }
    }

    selection.aggregate[[setting]] = list()
    ss = list()    
    for (v in vars) {
        ss[[v]] = apply(selection[[v]], 1, sum) / ncol(selection[[v]])
        selection.aggregate[[setting]][[v]] = ss[[v]]
    }
}


for (i in settings) {
    pdf(paste("figures/simulation/28-", i, "-profile-selection.pdf", sep=''))
    #dev.new()
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=1, xlab="Y-location", ylab="selection frequency")
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X2']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=2, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X3']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=3, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X4']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=4, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X5']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=5, ann=F, xaxt='n', yaxt='n')
    legend(legend=c("X1", "X2", "X3", "X4", "X5"), lty=1:5, x="topleft", bty='n')
    #abline(h=0.95, lty=3, col='red')
    title(paste("Simulation setting ", i, sep=""))
    dev.off()
}


coverage.oracular.bootstrap[[8]][['X1']] = coverage.oracular.bootstrap[[8]][['X1']][,-70]
for (setting in 1:18) {
    pdf(paste("figures/X1-", cluster, "-", setting, ".pdf", sep=""), width=12, height=4)
    layout(matrix(1:3,1,3))
    gwr.matplot(matrix(apply(coverage.unshrunk.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,1), axes=FALSE, x.crit=0.95, cs1.crit=1, cs2.crit=1, cs3.crit=1)
    axis(1, at=c(0,30), labels=c(0,1))
    axis(2, at=c(0,30), labels=c(0,1))
    title("GAL - unshrunk bootstrap")
    gwr.matplot(matrix(apply(coverage.oracular.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,1), axes=FALSE, x.crit=0.95, cs1.crit=1, cs2.crit=1, cs3.crit=1)
    axis(1, at=c(0,30), labels=c(0,1))   
    axis(2, at=c(0,30), labels=c(0,1))
    title("Oracular - bootstrap")
    gwr.matplot(matrix(apply(coverage.unshrunk.bootstrap[[setting]][['X1']],1,mean)/apply(coverage.oracular.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,5), axes=FALSE, x.crit=1, cs1.crit=1, cs2.crit=1, cs3.crit=1)
    axis(1, at=c(0,30), labels=c(0,1))   
    axis(2, at=c(0,30), labels=c(0,1))
    title("relative efficiency")
    dev.off()    
}


for (setting in 1:18) {
    pdf(paste("figures/selection-X1-", cluster, "-", setting, ".pdf", sep=""), width=12, height=4)
    layout(matrix(1:3,1,3))
    gwr.matplot(matrix(apply(coverage.unshrunk.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,1), axes=FALSE, x.crit=0.95, cs1.crit=1, cs2.crit=1, cs3.crit=1)
    axis(1, at=c(0,30), labels=c(0,1))
    axis(2, at=c(0,30), labels=c(0,1))
    title("GAL - unshrunk bootstrap")
    gwr.matplot(matrix(apply(coverage.oracular.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,1), axes=FALSE, x.crit=0.95, cs1.crit=1, cs2.crit=1, cs3.crit=1)
    axis(1, at=c(0,30), labels=c(0,1))   
    axis(2, at=c(0,30), labels=c(0,1))
    title("Oracular - bootstrap")
    gwr.matplot(matrix(apply(coverage.unshrunk.bootstrap[[setting]][['X1']],1,mean)/apply(coverage.oracular.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,5), axes=FALSE, x.crit=1, cs1.crit=1, cs2.crit=1, cs3.crit=1)
    axis(1, at=c(0,30), labels=c(0,1))   
    axis(2, at=c(0,30), labels=c(0,1))
    title("relative efficiency")
    dev.off()    
}


    gwr.matplot(matrix(apply(coverage.unshrunk.bootstrap[[setting]][['X1']],1,mean)/apply(coverage.oracular.bootstrap[[setting]][['X1']],1,mean),30,30), ylab="", xlab="", yrev=FALSE, cs1 = c(1, 0), cs2 = c(0.55, 0), cs3 = c(0, 0.55), border=NA, show.legend=TRUE, xrange=c(0,5), axes=FALSE, x.crit=1, cs1.crit=1, cs2.crit=1, cs3.crit=1)
