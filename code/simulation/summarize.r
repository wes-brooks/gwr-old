library(plotrix)

vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
params = c('bw', 'sigma2', 'loss.local', 's')

source('code/matplot.r')

#args = commandArgs(trailingOnly=TRUE)
#cluster = as.integer(args[1])
#cluster = 'NA'
cluster = 28

B = 100
N = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(c(0.03, 0.1), each=9)
rho = rep(rep(c(0, 0.5, 0.8), each=3), times=2)
sigma.tau = rep(c(0, 0.03, 0.1), times=6)
#function.type = rep(c("step", "gradient"), times=18)
#params = data.frame(tau, rho, sigma.tau, function.type)
params = data.frame(tau, rho, sigma.tau)

N = 30
B = list()
settings = 1:18

coord = seq(0, 1, length.out=N)
B[['(Intercept)']] = rep(0, N**2)
B[['X1']] = as.vector(matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N))
B[['X2']] = rep(0, N**2)
B[['X3']] = rep(0, N**2)
B[['X4']] = rep(0, N**2)
B[['X5']] = rep(0, N**2)

mean.coverage.ub = data.frame()
mean.coverage.b = data.frame()
mean.coverage.se = data.frame()
mean.selection = data.frame()

mean.coverage.oracular.b = data.frame()
mean.coverage.oracular.se = data.frame()


coverage.ub.aggregate = list()
coverage.b.aggregate = list()
coverage.se.aggregate = list()
selection.aggregate = list()

coverage.oracular.b.aggregate = list()
coverage.oracular.se.aggregate = list()

coverage.bootstrap = list()
coverage.unshrunk.bootstrap = list()
coverage.se = list()
selection = list()

coverage.oracular.bootstrap = list()
coverage.oracular.se = list()

for (setting in settings) {
    coverage.bootstrap[[setting]] = list()
    coverage.unshrunk.bootstrap[[setting]] = list()
    coverage.se[[setting]] = list()
    selection[[setting]] = list()

    coverage.oracular.bootstrap[[setting]] = list()
    coverage.oracular.se[[setting]] = list()

    #vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
    vars = c("X1")

    nsims = ifelse(setting==18,99,100)
    for (k in 1:nsims) {
        sim = (setting-1)*100 + k

        #Import our coefficient estimates
        filename = paste("output/output/CoefEstimates.", cluster, ".", sim, ".csv", sep="")
        filenameUnshrunk = paste("output/output/CoefEstimatesUnshrunk.", cluster, ".", sim, ".csv", sep="")
        filenameSEs = paste("output/output/CoefSEsUnshrunk.", cluster, ".", sim, ".csv", sep="")

        filenameCoefsOracular = paste("output/output/CoefEstimatesOracular.", cluster, ".", sim, ".csv", sep="")
        filenameSEsOracular = paste("output/output/CoefSEsOracular.", cluster, ".", sim, ".csv", sep="")

        estimates = read.csv(filename, header=TRUE)
        estimates.unshrunk = read.csv(filenameUnshrunk, header=TRUE)
        estimates.se = read.csv(filenameSEs, header=TRUE)

        estimates.oracular = read.csv(filenameCoefsOracular, header=TRUE)
        estimates.se.oracular = read.csv(filenameSEsOracular, header=TRUE)

        colnames(estimates) = vars
        colnames(estimates.unshrunk) = vars
        colnames(estimates.se) = vars

        colnames(estimates.oracular) = vars
        colnames(estimates.se.oracular) = vars
    
        for (v in vars) {
            #Calculate the coverage of each 95% CI 
            filename = paste("output/output/", v, ".", cluster, ".", sim, ".bootstrap.csv", sep="")
            bootstraps = as.matrix(read.csv(filename, header=TRUE))
            
            filename2 = paste("output/output/", v, ".", cluster, ".", sim, ".unshrunk-bootstrap.csv", sep="")
            unshrunk.bootstraps = as.matrix(read.csv(filename2, header=TRUE))

            filename = paste("output/output/", v, ".", cluster, ".", sim, ".OracularBootstrap.csv", sep="")
            oracularBootstraps = as.matrix(read.csv(filename, header=TRUE))
            
            CI.ub = t(apply(unshrunk.bootstraps, 1, function(x) {sort(x)[c(4, 98)]}))
            CI.b = t(apply(bootstraps, 1, function(x) {sort(x)[c(4, 98)]}))
            CI.se = cbind(estimates.unshrunk[,v]-1.96*estimates.se[,v], estimates.unshrunk[,v]+1.96*estimates.se[,v])

            CI.oracular.b = t(apply(oracularBootstraps, 1, function(x) {sort(x)[c(4, 98)]}))
            CI.oracular.se = cbind(estimates.oracular[,v]-1.96*estimates.se.oracular[,v], estimates.oracular[,v]+1.96*estimates.se.oracular[,v])
            
            if (k==1) {
                coverage.bootstrap[[setting]][[v]] = as.matrix(ifelse(B[[v]] < CI.b[,1] | B[[v]] > CI.b[,2],0,1))
            } else {
                coverage.bootstrap[[setting]][[v]] = cbind(coverage.bootstrap[[setting]][[v]], as.matrix(ifelse(B[[v]] < CI.b[,1] | B[[v]] > CI.b[,2],0,1)))
            }

            if (k==1) {
                coverage.unshrunk.bootstrap[[setting]][[v]] = as.matrix(ifelse(B[[v]] < CI.ub[,1] | B[[v]] > CI.ub[,2],0,1))
            } else {
                coverage.unshrunk.bootstrap[[setting]][[v]] = cbind(coverage.unshrunk.bootstrap[[setting]][[v]], as.matrix(ifelse(B[[v]] < CI.ub[,1] | B[[v]] > CI.ub[,2],0,1)))
            }
    
            if (k==1) {
                coverage.se[[setting]][[v]] = as.matrix(ifelse(B[[v]] < CI.se[,1] | B[[v]] > CI.se[,2],0,1))
            } else {
                coverage.se[[setting]][[v]] = cbind(coverage.se[[setting]][[v]], as.matrix(ifelse(B[[v]] < CI.se[,1] | B[[v]] > CI.se[,2],0,1)))
            }



            if (k==1) {
                coverage.oracular.bootstrap[[setting]][[v]] = as.matrix(ifelse(B[[v]] < CI.oracular.b[,1] | B[[v]] > CI.oracular.b[,2],0,1))
            } else {
                coverage.oracular.bootstrap[[setting]][[v]] = cbind(coverage.oracular.bootstrap[[setting]][[v]], as.matrix(ifelse(B[[v]] < CI.oracular.b[,1] | B[[v]] > CI.oracular.b[,2],0,1)))
            }
    
            if (k==1) {
                coverage.oracular.se[[setting]][[v]] = as.matrix(ifelse(B[[v]] < CI.oracular.se[,1] | B[[v]] > CI.oracular.se[,2],0,1))
            } else {
                coverage.oracular.se[[setting]][[v]] = cbind(coverage.oracular.se[[setting]][[v]], as.matrix(ifelse(B[[v]] < CI.oracular.se[,1] | B[[v]] > CI.oracular.se[,2],0,1)))
            }

            #Calculate how often each variable was selected for inclusion in the local models
            col = which(colnames(estimates) == v)
            if (k==1) {
                selection[[setting]][[v]] = as.matrix(ifelse(estimates[,col]==0, 0, 1))
            } else {
                selection[[setting]][[v]] = cbind(selection[[setting]][[v]], as.matrix(ifelse(estimates[,col]==0, 0, 1)))
            }
        }
    }

    cb = list()
    cub = list()
    cs = list()
    ss = list()

    cbo = list()
    cso = list()

    coverage.ub.aggregate[[setting]] = list()
    coverage.b.aggregate[[setting]] = list()
    coverage.se.aggregate[[setting]] = list()
    selection.aggregate[[setting]] = list()

    coverage.oracular.b.aggregate[[setting]] = list()
    coverage.oracular.se.aggregate[[setting]] = list()

    for (v in c("X1")) { #vars) {
        cb[[v]] = apply(coverage.bootstrap[[setting]][[v]], 1, sum) / ncol(coverage.bootstrap[[setting]][[v]])
        cub[[v]] = apply(coverage.unshrunk.bootstrap[[setting]][[v]], 1, sum) / ncol(coverage.unshrunk.bootstrap[[setting]][[v]])
        cs[[v]] = apply(coverage.se[[setting]][[v]], 1, sum) / ncol(coverage.se[[setting]][[v]])
        ss[[v]] = apply(selection[[setting]][[v]], 1, sum) / ncol(selection[[setting]][[v]])

        cbo[[v]] = apply(coverage.oracular.bootstrap[[setting]][[v]], 1, sum) / ncol(coverage.oracular.bootstrap[[setting]][[v]])
        cso[[v]] = apply(coverage.oracular.se[[setting]][[v]], 1, sum) / ncol(coverage.oracular.se[[setting]][[v]])       

#        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".unshrunk_bootstrap_coverage.pdf", sep=""))
#        gwr.matplot(matrix(cub[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#        dev.off()

#        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".oracular_bootstrap_coverage.pdf", sep=""))
#        gwr.matplot(matrix(cbo[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#        dev.off()

#        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".oracular_se_coverage.pdf", sep=""))
#        gwr.matplot(matrix(cso[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#        dev.off()
        
#         pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".bootstrap_coverage.pdf", sep=""))
#         gwr.matplot(matrix(cb[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#         #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#         dev.off()
#     
#        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".se_coverage.pdf", sep=""))
#        gwr.matplot(matrix(cs[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#        dev.off()
# 
#         pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".selection.pdf", sep=""))
#         gwr.matplot(matrix(ss[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#         #title(main=paste("Frequency that ", v, "is selected", sep=""))
#         dev.off()

        coverage.ub.aggregate[[setting]][[v]] = cub[[v]]
        coverage.b.aggregate[[setting]][[v]] = cb[[v]]
        coverage.se.aggregate[[setting]][[v]] = cs[[v]]
        selection.aggregate[[setting]][[v]] = ss[[v]]
    
        coverage.oracular.b.aggregate[[setting]][[v]] = cbo[[v]]
        coverage.oracular.se.aggregate[[setting]][[v]] = cso[[v]]
    }



#     vars = c("X1", "X2", "X3", "X4", "X5")
#     xx = c(0, 30)
#     yy = c(0,1)
#     
#     pdf(paste("figures/simulation/", cluster, ".", setting, ".profile_selection.pdf", sep=""))
#     plot(apply(matrix(ss[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="selection frequency")
#     for (k in 2:length(vars)) {
#         v = vars[k]
#         par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
#         plot(apply(matrix(ss[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
#     }
#     legend(vars, lty=1:length(vars), bty='n', x='topleft')
#     dev.off()
#     
#     
#     pdf(paste("figures/simulation/", cluster, ".", setting, ".profile_bootstrap_coverage.pdf", sep=""))
#     vars = c("X1", "X2", "X3", "X4", "X5", "(Intercept)")
#     plot(apply(matrix(cs[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
#     for (k in 2:length(vars)) {
#         v = vars[k]
#         par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
#         plot(apply(matrix(cs[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
#     }
#     abline(h=0.95, lty=3, col='red')
#     legend(vars, lty=1:length(vars), bty='n', x='bottomleft')
#     dev.off()
#     
#     
#     pdf(paste("figures/simulation/", cluster, ".", setting, ".profile_se_coverage.pdf", sep=""))
#     plot(apply(matrix(cb[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
#     for (k in 2:length(vars)) {
#         v = vars[k]
#         par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
#         plot(apply(matrix(cb[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
#     }
#     abline(h=0.95, lty=3, col='red')
#     legend(vars, lty=1:length(vars), bty='n', x='bottomleft')
#     dev.off()
#   

#    cbv = vector()
#    cubv = vector()
#    csv = vector()
#    sv = vector()
#
#    cbov = vector()
#    csov = vector()
#
#    for (k in 1:length(vars)) {
#        v = vars[k]
#        csv = c(csv, mean(cs[[v]]))
#        cbv = c(cbv, mean(cb[[v]]))
#        cubv = c(cubv, mean(cub[[v]]))
#
#        csov = c(csov, mean(cso[[v]]))
#        cbov = c(cbov, mean(cbo[[v]]))
#        
#        if (v!="(Intercept)") {
#            sv = c(sv, mean(ifelse(B[[v]]==0, 1-ss[[v]], ss[[v]])))
#        }
#    }
#    
#    mean.coverage.ub = rbind(mean.coverage.ub, cubv)
#    mean.coverage.b = rbind(mean.coverage.b, cbv)
#    mean.coverage.se = rbind(mean.coverage.se, csv)
#    mean.selection = rbind(mean.selection, sv)
#
#    mean.coverage.oracular.b = rbind(mean.coverage.oracular.b, cbov)
#    mean.coverage.oracular.se = rbind(mean.coverage.oracular.se, csov)



}

for (i in settings) {
    pdf(paste("figures/simulation/28-", i, "-profile-coverage.pdf", sep=''))
    #dev.new()
    plot(x=seq(0,1,length.out=30), y=apply(matrix(coverage.b.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=1, xlab="Y-location", ylab="coverage")
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(coverage.ub.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=2, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(coverage.se.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=3, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(coverage.oracular.b.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=4, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(coverage.oracular.se.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=5, ann=F, xaxt='n', yaxt='n')
    legend(legend=c("LASSO/bootstrap", "LASSO/unshrunk bootstrap", "LASSO/SE", "Oracular/bootstrap", "Oracular/SE"), lty=1:5, x="bottomright", bty='n')
    abline(h=0.95, lty=3, col='red')
    title(paste("Simulation setting ", i, sep=""))
    dev.off()
}

for (i in settings) {
    #pdf(paste("figures/simulation/28-", i, "-profile-selection.pdf", sep=''))
    dev.new()
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X1']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=1, xlab="Y-location", ylab="selection frequency")
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X2']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=2, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X3']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=3, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X4']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=4, ann=F, xaxt='n', yaxt='n')
    par(new=TRUE)
    plot(x=seq(0,1,length.out=30), y=apply(matrix(selection.aggregate[[i]][['X5']],ncol=30,nrow=30), 1, mean), type='l', xlim=c(0,1), ylim=c(0,1), bty='n', lty=5, ann=F, xaxt='n', yaxt='n')
    legend(legend=c("X1", "X2", "LASSO/SE", "X3", "X4"), lty=1:5, x="topleft", bty='n')
    #abline(h=0.95, lty=3, col='red')
    title(paste("Simulation setting ", i, sep=""))
    #dev.off()
}

#    pdf(paste("figures/simulation/", cluster, ".profile_unshrunk_bootstrap_coverage.pdf", sep=""))
#    vars = c("X1", "X2", "X3", "X4", "X5", "(Intercept)")
#    plot(apply(matrix(cub[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
#    for (setting in settings) {
#        v = vars[k]
#        par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
#        plot(apply(matrix(cub[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
#    }
#    abline(h=0.95, lty=3, col='red')
#    legend(vars, lty=1:length(vars), bty='n', x='bottomleft')
#    dev.off()
#    
#
#    vars = c("X1", "X2", "X3", "X4", "X5", "(Intercept)")
#
#for (v in vars) {
#    pdf(paste("figures/simulation/", cluster, ".", v, ".profile_oracular_bootstrap_coverage.pdf", sep=""))
#    plot(apply(matrix(cbo[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
#    for (setting in settings) {
#        v = vars[k]
#        par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
#        plot(apply(matrix(cbo[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
#    }
#    abline(h=0.95, lty=3, col='red')
#    legend(vars, lty=1:length(vars), bty='n', x='bottomleft')
#    dev.off()
#}
#
#
#
#colnames(mean.coverage.ub) = vars
#colnames(mean.coverage.b) = vars
#colnames(mean.coverage.se) = vars
#colnames(mean.selection) = vars[2:6]
#colnames(mean.coverage.oracular.b) = vars
#colnames(mean.coverage.oracular.se) = vars

t.x = list()
#t.x[[1]] = 1:16
#t.x[[2]] = 17:36
#t.x[[1]] = 1:9
#t.x[[2]] = 10:18
t.x[[1]] = 1:4
#t.x[[2]] = 10:18

rho = list()
#rho[[1]] = c(1:6, 19:24)
#rho[[2]] = c(7:12, 25:30)
#rho[[3]] = c(13:18, 31:36)
#rho[[1]] = c(1:6)
#rho[[2]] = c(7:12)
#rho[[3]] = c(13:18)
rho[[1]] = c(1:4)

t.e = list()
#t.e[[1]] = c(1:2, 7:8, 13:14, 19:20, 25:26, 31:32)
#t.e[[2]] = c(3:4, 9:10, 15:16, 21:22, 27:28, 33:34)
#t.e[[3]] = c(5:6, 11:12, 17:18, 23:24, 29:30, 35:36)
#t.e[[1]] = c(1:2, 7:8, 13:14)
#t.e[[2]] = c(3:4, 9:10, 15:16)
#t.e[[3]] = c(5:6, 11:12, 17:18)
t.e[[1]] = c(1:2)
t.e[[2]] = c(3:4)

function.type = list()
function.type[['step']] = (1:18) - 1
function.type[['gradient']] = (1:18)*2


#contrasts = c("t.x", "rho", "t.e", "function.type")
contrasts = c("t.x", "rho", "t.e")
for (contrast in contrasts) {
    contr = get(contrast)
    dev.new()
    #pdf(paste("figures/simulation/", cluster, ".", contrast, ".profile_bootstrap_coverage.pdf", sep=""))
    for (i in 1:length(contr)) {
        element = contr[[i]]

        #Isolate data for the boxplot
        plotdata = data.frame()
        for (j in element) {
            plotdata = rbind(plotdata, coverage.ub.aggregate[[j]][['X1']])
        }

        boxplot(plotdata, col=i)
    }

        plot(apply(matrix(cs[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")


        for (k in 2:length(vars)) {
            v = vars[k]
            par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
            plot(apply(matrix(cs[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
        }
        abline(h=0.95, lty=3, col='red')
        legend(vars, lty=1:length(vars), bty='n', x='bottomleft')
        dev.off()
    }
}

t.x.selection.table = rbind(apply(mean.selection[t.x[[1]],],2,mean), apply(mean.selection[t.x[[2]],],2,mean))
t.x.b.coverage.table = rbind(apply(mean.coverage.b[t.x[[1]],],2,mean), apply(mean.coverage.b[t.x[[2]],],2,mean))
t.x.ub.coverage.table = rbind(apply(mean.coverage.ub[t.x[[1]],],2,mean), apply(mean.coverage.ub[t.x[[2]],],2,mean))
t.x.s.coverage.table = rbind(apply(mean.coverage.se[t.x[[1]],],2,mean), apply(mean.coverage.se[t.x[[2]],],2,mean))
t.x.o.b.coverage.table = rbind(apply(mean.coverage.oracular.b[t.x[[1]],],2,mean), apply(mean.coverage.oracular.b[t.x[[2]],],2,mean))
t.x.o.se.coverage.table = rbind(apply(mean.coverage.oracular.se[t.x[[1]],],2,mean), apply(mean.coverage.oracular.se[t.x[[2]],],2,mean))
rownames(t.x.selection.table) = rownames(t.x.o.b.coverage.table) = rownames(t.x.o.se.coverage.table) = rownames(t.x.b.coverage.table) = rownames(t.x.ub.coverage.table) = rownames(t.x.s.coverage.table) = c("0.1", "0.8")
xtable(t.x.selection.table)
xtable(t.x.b.coverage.table)
xtable(t.x.ub.coverage.table)
xtable(t.x.s.coverage.table)
xtable(t.x.o.b.coverage.table)
xtable(t.x.o.se.coverage.table)

t.e.selection.table = rbind(apply(mean.selection[t.e[[1]],],2,mean), apply(mean.selection[t.e[[2]],],2,mean), apply(mean.selection[t.e[[3]],],2,mean))
t.e.b.coverage.table = rbind(apply(mean.coverage.b[t.e[[1]],],2,mean), apply(mean.coverage.b[t.e[[2]],],2,mean), apply(mean.coverage.b[t.e[[3]],],2,mean))
t.e.ub.coverage.table = rbind(apply(mean.coverage.ub[t.e[[1]],],2,mean), apply(mean.coverage.ub[t.e[[2]],],2,mean), apply(mean.coverage.ub[t.e[[3]],],2,mean))
t.e.s.coverage.table = rbind(apply(mean.coverage.se[t.e[[1]],],2,mean), apply(mean.coverage.se[t.e[[2]],],2,mean), apply(mean.coverage.se[t.e[[3]],],2,mean))
t.e.o.b.coverage.table = rbind(apply(mean.coverage.oracular.b[t.e[[1]],],2,mean), apply(mean.coverage.oracular.b[t.e[[2]],],2,mean), apply(mean.coverage.oracular.b[t.e[[3]],],2,mean))
t.e.o.se.coverage.table = rbind(apply(mean.coverage.oracular.se[t.e[[1]],],2,mean), apply(mean.coverage.oracular.se[t.e[[2]],],2,mean), apply(mean.coverage.oracular.se[t.e[[3]],],2,mean))
rownames(t.e.selection.table) = rownames(t.e.o.b.coverage.table) = rownames(t.e.o.se.coverage.table) = rownames(t.e.b.coverage.table) = rownames(t.e.ub.coverage.table) = rownames(t.e.s.coverage.table) = c("0", "0.1", "0.8")
xtable(t.e.selection.table)
xtable(t.e.b.coverage.table)
xtable(t.e.ub.coverage.table)
xtable(t.e.s.coverage.table)
xtable(t.e.o.b.coverage.table)
xtable(t.e.o.se.coverage.table)

rho.selection.table = rbind(apply(mean.selection[rho[[1]],],2,mean), apply(mean.selection[rho[[2]],],2,mean), apply(mean.selection[rho[[3]],],2,mean))
rho.b.coverage.table = rbind(apply(mean.coverage.b[rho[[1]],],2,mean), apply(mean.coverage.b[rho[[2]],],2,mean), apply(mean.coverage.b[rho[[3]],],2,mean))
rho.ub.coverage.table = rbind(apply(mean.coverage.ub[rho[[1]],],2,mean), apply(mean.coverage.ub[rho[[2]],],2,mean), apply(mean.coverage.ub[rho[[3]],],2,mean))
rho.s.coverage.table = rbind(apply(mean.coverage.se[rho[[1]],],2,mean), apply(mean.coverage.se[rho[[2]],],2,mean), apply(mean.coverage.se[rho[[3]],],2,mean))
rho.o.b.coverage.table = rbind(apply(mean.coverage.oracular.b[rho[[1]],],2,mean), apply(mean.coverage.oracular.b[rho[[2]],],2,mean), apply(mean.coverage.oracular.b[rho[[3]],],2,mean))
rho.o.se.coverage.table = rbind(apply(mean.coverage.oracular.se[rho[[1]],],2,mean), apply(mean.coverage.oracular.se[rho[[2]],],2,mean), apply(mean.coverage.oracular.se[rho[[3]],],2,mean))
rownames(rho.selection.table) = rownames(rho.o.b.coverage.table) = rownames(rho.o.se.coverage.table) = rownames(rho.b.coverage.table) = rownames(rho.ub.coverage.table) = rownames(rho.s.coverage.table) = c("0", "0.5", "0.8")
xtable(rho.selection.table)
xtable(rho.b.coverage.table)
xtable(rho.ub.coverage.table)
xtable(rho.s.coverage.table)
xtable(rho.o.b.coverage.table)
xtable(rho.o.se.coverage.table)

func.selection.table = rbind(apply(mean.selection[function.type[['gradient']],],2,mean), apply(mean.selection[function.type[['step']],],2,mean))
func.b.coverage.table = rbind(apply(mean.coverage.b[function.type[['gradient']],],2,mean), apply(mean.coverage.b[function.type[['step']],],2,mean))
func.ub.coverage.table = rbind(apply(mean.coverage.ub[function.type[['gradient']],],2,mean), apply(mean.coverage.ub[function.type[['step']],],2,mean))
func.s.coverage.table = rbind(apply(mean.coverage.se[function.type[['gradient']],],2,mean), apply(mean.coverage.se[function.type[['step']],],2,mean))
func.o.b.coverage.table = rbind(apply(mean.coverage.oracular.b[function.type[['gradient']],],2,mean), apply(mean.coverage.oracular.b[function.type[['step']],],2,mean))
func.o.se.coverage.table = rbind(apply(mean.coverage.oracular.se[function.type[['gradient']],],2,mean), apply(mean.coverage.oracular.se[function.type[['step']],],2,mean))
rownames(func.selection.table) = rownames(func.o.b.coverage.table) = rownames(func.o.se.coverage.table) = rownames(func.b.coverage.table) = rownames(func.ub.coverage.table) = rownames(func.s.coverage.table) = c("gradient", "step")
xtable(func.selection.table)
xtable(func.b.coverage.table)
xtable(func.ub.coverage.table)
xtable(func.s.coverage.table)
xtable(func.o.b.coverage.table)
xtable(func.o.se.coverage.table)