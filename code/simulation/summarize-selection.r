library(plotrix)

vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
params = c('bw', 'sigma2', 'loss.local', 's')

source('code/matplot.r')

#args = commandArgs(trailingOnly=TRUE)
#cluster = as.integer(args[1])
#cluster = 'NA'
cluster = 43

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


selection.aggregate = list()


selection = list()

selection.precon = list()

for (setting in settings) {
	cat(paste("Begin setting ", setting, ".\n", sep=""))
    selection[[setting]] = list()
    

	selection.precon[[setting]] = list()
    vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
    vv = c("X1")

    nsims = 100
    for (k in 1:nsims) {
        sim = (setting-1)*100 + k - 1
        cat(paste("Begin simulation ", sim, ".\n", sep=""))      
        
        #Import our coefficient estimates
        filename = paste("output/CoefEstimates.", cluster, ".", sim, ".csv", sep="")
        filename.precon = paste("output/CoefEstimates.", cluster, ".", sim, ".precon.csv", sep="")

        estimates = read.csv(filename, header=TRUE)
		estimates.precon = read.csv(filename.precon, header=TRUE)

        colnames(estimates) = vars
        colnames(estimates.precon) = vars
        
        for (v in vars) {
        	cat(paste("Begin variable ", v, "\n", sep=""))

            #Calculate how often each variable was selected for inclusion in the local models
            col = which(colnames(estimates) == v)
            if (k==1) {
                selection[[setting]][[v]] = as.matrix(ifelse(estimates[,col]==0, 0, 1))
            } else {
                selection[[setting]][[v]] = cbind(selection[[setting]][[v]], as.matrix(ifelse(estimates[,col]==0, 0, 1)))
            }
            
            #Calculate how often each variable was selected for inclusion in the local models
            col = which(colnames(estimates.precon) == v)
            if (k==1) {
                selection.precon[[setting]][[v]] = as.matrix(ifelse(estimates.precon[,col]==0, 0, 1))
            } else {
                selection.precon[[setting]][[v]] = cbind(selection.precon[[setting]][[v]], as.matrix(ifelse(estimates.precon[,col]==0, 0, 1)))
            }
        }
    }
}

	cat(paste("Now aggregate.\n", sep=""))
    cb = list()
    cub = list()
    cs = list()
    ss = list()
    
	cbp = list()
    cubp = list()
    csp = list()
    ssp = list()

    cbo = list()
    cso = list()

    coverage.ub.precon.aggregate[[setting]] = list()
    coverage.b.precon.aggregate[[setting]] = list()
    coverage.se.precon.aggregate[[setting]] = list()
    selection.precon.aggregate[[setting]] = list()

    coverage.ub.aggregate[[setting]] = list()
    coverage.b.aggregate[[setting]] = list()
    coverage.se.aggregate[[setting]] = list()
    selection.aggregate[[setting]] = list()

    coverage.oracular.b.aggregate[[setting]] = list()
    coverage.oracular.se.aggregate[[setting]] = list()

    for (v in c("X1")) { #vars) {
    	cat(paste("Aggregating ", v, ".\n", sep=""))
        cb[[v]] = apply(coverage.bootstrap[[setting]][[v]], 1, sum) / ncol(coverage.bootstrap[[setting]][[v]])
        cub[[v]] = apply(coverage.unshrunk.bootstrap[[setting]][[v]], 1, sum) / ncol(coverage.unshrunk.bootstrap[[setting]][[v]])
        cs[[v]] = apply(coverage.se[[setting]][[v]], 1, sum) / ncol(coverage.se[[setting]][[v]])
        ss[[v]] = apply(selection[[setting]][[v]], 1, sum) / ncol(selection[[setting]][[v]])
        
        cbp[[v]] = apply(coverage.bootstrap.precon[[setting]][[v]], 1, sum) / ncol(coverage.bootstrap.precon[[setting]][[v]])
        cubp[[v]] = apply(coverage.unshrunk.bootstrap.precon[[setting]][[v]], 1, sum) / ncol(coverage.unshrunk.bootstrap.precon[[setting]][[v]])
        csp[[v]] = apply(coverage.se.precon[[setting]][[v]], 1, sum) / ncol(coverage.se.precon[[setting]][[v]])
        ssp[[v]] = apply(selection.precon[[setting]][[v]], 1, sum) / ncol(selection.precon[[setting]][[v]])

        cbo[[v]] = apply(coverage.oracular.bootstrap[[setting]][[v]], 1, sum) / ncol(coverage.oracular.bootstrap[[setting]][[v]])
        cso[[v]] = apply(coverage.oracular.se[[setting]][[v]], 1, sum) / ncol(coverage.oracular.se[[setting]][[v]])      

        coverage.ub.aggregate[[setting]][[v]] = cub[[v]]
        coverage.b.aggregate[[setting]][[v]] = cb[[v]]
        coverage.se.aggregate[[setting]][[v]] = cs[[v]]
        selection.aggregate[[setting]][[v]] = ss[[v]]
        
		coverage.ub.precon.aggregate[[setting]][[v]] = cubp[[v]]
        coverage.b.precon.aggregate[[setting]][[v]] = cbp[[v]]
        coverage.se.precon.aggregate[[setting]][[v]] = csp[[v]]
        selection.precon.aggregate[[setting]][[v]] = ssp[[v]]
    
        coverage.oracular.b.aggregate[[setting]][[v]] = cbo[[v]]
        coverage.oracular.se.aggregate[[setting]][[v]] = cso[[v]]
    }
}