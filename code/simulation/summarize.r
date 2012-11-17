library(plotrix)

vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
params = c('bw', 'sigma2', 'loss.local', 's', 'sum.weights')

source('code/matplot.r')

args = commandArgs(trailingOnly=TRUE)
#cluster = as.integer(args[1])
cluster = 15

B = 100
N = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(c(0.1, 0.8), each=18)
rho = rep(rep(c(0, 0.5, 0.8), each=6), times=2)
sigma.tau = rep(rep(c(0, 0.1, 0.8), each=2), times=6)
function.type = rep(c("step", "gradient"), times=18)
params = data.frame(tau, rho, sigma.tau, function.type)


N = 30
B = list()

coord = seq(0, 1, length.out=N)
B[['(Intercept)']] = rep(0, N**2)
B[['X2']] = rep(0, N**2)
B[['X3']] = rep(0, N**2)
B[['X4']] = rep(0, N**2)
B[['X5']] = rep(0, N**2)

mean.coverage.b = data.frame()
mean.coverage.s = data.frame()
mean.selection = data.frame()

for (setting in 1:36) {
    coverage.bootstrap = list()
    coverage.se = list()
    selection = list()

    vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')

    if (setting %% 2 == 1) {
        B[['X1']] = as.vector(matrix(rep(ifelse(coord<=0.5, 0, 1), N), N, N))
    } else {
        B[['X1']] = as.vector(matrix(rep(1-coord, N), N, N))
    }

    nsims = ifelse(setting<36,100,99)
    for (k in 1:nsims) {
        sim = (setting-1)*100 + k

        #Import our coefficient estimates
        filename = paste("output/output/CoefEstimates.", cluster, ".", sim, ".csv", sep="")
        filenameUnshrunk = paste("output/output/CoefEstimatesUnshrunk.", cluster, ".", sim, ".csv", sep="")
        filenameSEs = paste("output/output/CoefSEsUnshrunk.", cluster, ".", sim, ".csv", sep="")

        estimates = read.csv(filename, header=TRUE)
        estimates.unshrunk = read.csv(filenameUnshrunk, header=TRUE)
        estimates.se = read.csv(filenameSEs, header=TRUE)

        colnames(estimates) = vars
        colnames(estimates.unshrunk) = vars
        colnames(estimates.se) = vars
    
        for (v in vars) {
            #Calculate the coverage of each 95% CI 
            filename = paste("output/output/", v, ".", cluster, ".", sim, ".bootstrap.csv", sep="")
            bootstraps = as.matrix(read.csv(filename, header=TRUE))
            
            filename2 = paste("output/output/", v, ".", cluster, ".", sim, ".unshrunk-bootstrap.csv", sep="")
            unshrunk.bootstraps = as.matrix(read.csv(filename2, header=TRUE))
            
            CI.ub = t(apply(unshrunk.bootstraps, 1, function(x) {sort(x)[c(4, 98)]}))
            CI.b = t(apply(bootstraps, 1, function(x) {sort(x)[c(4, 98)]}))

            CI.se = cbind(estimates.unshrunk[,v]-1.96*estimates.se[,v], estimates.unshrunk[,v]+1.96*estimates.se[,v])
            
            if (k==0) {
                coverage.bootstrap[[v]] = as.matrix(ifelse(B[[v]] < CI.b[,1] | B[[v]] > CI.b[,2],0,1))
            } else {
                coverage.bootstrap[[v]] = cbind(coverage.bootstrap[[v]], as.matrix(ifelse(B[[v]] < CI.b[,1] | B[[v]] > CI.b[,2],0,1)))
            }

            if (k==0) {
                coverage.unshrunk.bootstrap[[v]] = as.matrix(ifelse(B[[v]] < CI.ub[,1] | B[[v]] > CI.ub[,2],0,1))
            } else {
                coverage.unshrunk.bootstrap[[v]] = cbind(coverage.unshrunk.bootstrap[[v]], as.matrix(ifelse(B[[v]] < CI.ub[,1] | B[[v]] > CI.ub[,2],0,1)))
            }
    
            if (k==0) {
                coverage.se[[v]] = as.matrix(ifelse(B[[v]] < CI.se[,1] | B[[v]] > CI.se[,2],0,1))
            } else {
                coverage.se[[v]] = cbind(coverage.se[[v]], as.matrix(ifelse(B[[v]] < CI.se[,1] | B[[v]] > CI.se[,2],0,1)))
            }

            #Calculate how often each variable was selected for inclusion in the local models
            col = which(colnames(estimates) == v)
            if (k==0) {
                selection[[v]] = as.matrix(ifelse(estimates[,col]==0, 0, 1))
            } else {
                selection[[v]] = cbind(selection[[v]], as.matrix(ifelse(estimates[,col]==0, 0, 1)))
            }        
        }
    }

    cb = list()
    cub = list()
    cs = list()
    ss = list()
    for (v in vars) {
        cb[[v]] = apply(coverage.bootstrap[[v]], 1, sum) / ncol(coverage.bootstrap[[v]])
        cub[[v]] = apply(coverage.bootstrap[[v]], 1, sum) / ncol(coverage.bootstrap[[v]])
        cs[[v]] = apply(coverage.se[[v]], 1, sum) / ncol(coverage.se[[v]])
        ss[[v]] = apply(selection[[v]], 1, sum) / ncol(selection[[v]])

        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".unshrunk_bootstrap_coverage.pdf", sep=""))
        gwr.matplot(matrix(cb[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
        dev.off()
        
#         pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".bootstrap_coverage.pdf", sep=""))
#         gwr.matplot(matrix(cb[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#         #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#         dev.off()
#     
#         pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".se_coverage.pdf", sep=""))
#         gwr.matplot(matrix(cs[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#         #title(main=paste("Coverage of 95% CI for ", v, sep=""))
#         dev.off()
# 
#         pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".selection.pdf", sep=""))
#         gwr.matplot(matrix(ss[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
#         #title(main=paste("Frequency that ", v, "is selected", sep=""))
#         dev.off()
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
    pdf(paste("figures/simulation/", cluster, ".", setting, ".profile_unshrunk_bootstrap_coverage.pdf", sep=""))
    vars = c("X1", "X2", "X3", "X4", "X5", "(Intercept)")
    plot(apply(matrix(cub[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
    for (k in 2:length(vars)) {
        v = vars[k]
        par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
        plot(apply(matrix(cub[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
    }
    abline(h=0.95, lty=3, col='red')
    legend(vars, lty=1:length(vars), bty='n', x='bottomleft')
    dev.off()
    
    cbv = vector()
    cubv = vector()
    csv = vector()
    sv = vector()
    for (k in 1:length(vars)) {
        v = vars[k]
        csv = c(csv, mean(cs[[v]]))
        cbv = c(cbv, mean(cb[[v]]))
        cubv = c(cubv, mean(cub[[v]]))
        
        if (v!="(Intercept)") {
            sv = c(sv, mean(ifelse(B[[v]]==0, 1-ss[[v]], ss[[v]])))
        }
    }
    
    mean.coverage.ub = rbind(mean.coverage.b, cbv)
    mean.coverage.b = rbind(mean.coverage.b, cbv)
    mean.coverage.s = rbind(mean.coverage.s, csv)
    mean.selection = rbind(mean.selection, sv)
}

colnames(mean.coverage.ub) = vars
colnames(mean.coverage.b) = vars
colnames(mean.coverage.s) = vars
colnames(mean.selection) = vars[2:6]


t.x.1 = 1:16
t.x.2 = 17:36

rho.1 = c(1:6, 19:24)
rho.2 = c(7:12, 25:30)
rho.3 = c(13:18, 31:36)

t.e.1 = c(1:2, 7:8, 13:14, 19:20, 25:26, 31:32)
t.e.2 = c(3:4, 9:10, 15:16, 21:22, 27:28, 33:34)
t.e.3 = c(5:6, 11:12, 17:18, 23:24, 29:30, 35:36)

step = (1:18)*2 - 1
gradient = (1:18)*2


t.x.selection.table = rbind(apply(mean.selection[t.x.1,],2,mean), apply(mean.selection[t.x.2,],2,mean))
t.x.b.coverage.table = rbind(apply(mean.coverage.b[t.x.1,],2,mean), apply(mean.coverage.b[t.x.2,],2,mean))
t.x.ub.coverage.table = rbind(apply(mean.coverage.ub[t.x.1,],2,mean), apply(mean.coverage.ub[t.x.2,],2,mean))
t.x.s.coverage.table = rbind(apply(mean.coverage.s[t.x.1,],2,mean), apply(mean.coverage.s[t.x.2,],2,mean))
rownames(t.x.selection.table) = rownames(t.x.b.coverage.table) = rownames(t.x.ub.coverage.table) = rownames(t.x.s.coverage.table) = c("0.1", "0.8")
xtable(t.x.selection.table)
xtable(t.x.b.coverage.table)
xtable(t.x.ub.coverage.table)
xtable(t.x.s.coverage.table)

t.e.selection.table = rbind(apply(mean.selection[t.e.1,],2,mean), apply(mean.selection[t.e.2,],2,mean), apply(mean.selection[t.e.3,],2,mean))
t.e.b.coverage.table = rbind(apply(mean.coverage.b[t.e.1,],2,mean), apply(mean.coverage.b[t.e.2,],2,mean), apply(mean.coverage.b[t.e.3,],2,mean))
t.e.ub.coverage.table = rbind(apply(mean.coverage.ub[t.e.1,],2,mean), apply(mean.coverage.ub[t.e.2,],2,mean), apply(mean.coverage.ub[t.e.3,],2,mean))
t.e.s.coverage.table = rbind(apply(mean.coverage.s[t.e.1,],2,mean), apply(mean.coverage.s[t.e.2,],2,mean), apply(mean.coverage.s[t.e.3,],2,mean))
rownames(t.e.selection.table) = rownames(t.e.b.coverage.table) = rownames(t.e.ub.coverage.table) = rownames(t.e.s.coverage.table) = c("0", "0.1", "0.8")
xtable(t.e.selection.table)
xtable(t.e.b.coverage.table)
xtable(t.e.ub.coverage.table)
xtable(t.e.s.coverage.table)

rho.selection.table = rbind(apply(mean.selection[rho.1,],2,mean), apply(mean.selection[rho.2,],2,mean), apply(mean.selection[rho.3,],2,mean))
rho.b.coverage.table = rbind(apply(mean.coverage.b[rho.1,],2,mean), apply(mean.coverage.b[rho.2,],2,mean), apply(mean.coverage.b[rho.3,],2,mean))
rho.ub.coverage.table = rbind(apply(mean.coverage.ub[rho.1,],2,mean), apply(mean.coverage.ub[rho.2,],2,mean), apply(mean.coverage.ub[rho.3,],2,mean))
rho.s.coverage.table = rbind(apply(mean.coverage.s[rho.1,],2,mean), apply(mean.coverage.s[rho.2,],2,mean), apply(mean.coverage.s[rho.3,],2,mean))
rownames(rho.selection.table) = rownames(rho.b.coverage.table) = rownames(rho.ub.coverage.table) = rownames(rho.s.coverage.table) = c("0", "0.5", "0.8")
xtable(rho.selection.table)
xtable(rho.b.coverage.table)
xtable(rho.ub.coverage.table)
xtable(rho.s.coverage.table)

func.selection.table = rbind(apply(mean.selection[gradient,],2,mean), apply(mean.selection[step,],2,mean))
func.b.coverage.table = rbind(apply(mean.coverage.b[gradient,],2,mean), apply(mean.coverage.b[step,],2,mean))
func.ub.coverage.table = rbind(apply(mean.coverage.ub[gradient,],2,mean), apply(mean.coverage.ub[step,],2,mean))
func.s.coverage.table = rbind(apply(mean.coverage.s[gradient,],2,mean), apply(mean.coverage.s[step,],2,mean))
rownames(func.selection.table) = rownames(func.b.coverage.table) = rownames(func.ub.coverage.table) = rownames(func.s.coverage.table) = c("gradient", "step")
xtable(func.selection.table)
xtable(func.b.coverage.table)
xtable(func.ub.coverage.table)
xtable(func.s.coverage.table)