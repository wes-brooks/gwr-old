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
B[['X1']] = as.vector(matrix(rep(ifelse(coord<=0.5, 0, 1), N), N, N))
B[['X2']] = rep(0, N**2)
B[['X3']] = rep(0, N**2)
B[['X4']] = rep(0, N**2)
B[['X5']] = rep(0, N**2)

coverage.bootstrap = list()
coverage.se = list()
selection = list()

for (setting in 1:1) {
    for (k in 1:100) {
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

            CI.b = t(apply(bootstraps, 1, function(x) {sort(x)[c(4, 98)]}))
            CI.se = cbind(estimates.unshrunk[,v]-1.96*estimates.se[,v], estimates.unshrunk[,v]+1.96*estimates.se[,v])
            
            if (k==0) {
                coverage.bootstrap[[v]] = as.matrix(ifelse(B[[v]] < CI.b[,1] | B[[v]] > CI.b[,2],0,1))
            } else {
                coverage.bootstrap[[v]] = cbind(coverage.bootstrap[[v]], as.matrix(ifelse(B[[v]] < CI.b[,1] | B[[v]] > CI.b[,2],0,1)))
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
    cs = list()
    ss = list()
    for (v in vars) {
        cb[[v]] = apply(coverage.bootstrap[[v]], 1, sum) / ncol(coverage.bootstrap[[v]])
        cs[[v]] = apply(coverage.se[[v]], 1, sum) / ncol(coverage.se[[v]])
        ss[[v]] = apply(selection[[v]], 1, sum) / ncol(selection[[v]])
        
        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".bootstrap_coverage.pdf", sep=""))
        gwr.matplot(matrix(cb[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
        dev.off()
    
        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".se_coverage.pdf", sep=""))
        gwr.matplot(matrix(cs[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
        #title(main=paste("Coverage of 95% CI for ", v, sep=""))
        dev.off()

        pdf(paste("figures/simulation/", v, ".", cluster, ".", setting, ".selection.pdf", sep=""))
        gwr.matplot(matrix(ss[[v]], N, N), c(0,1), c(0,1), c(0,1), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=c(0,1))
        #title(main=paste("Frequency that ", v, "is selected", sep=""))
        dev.off()
    }

}
