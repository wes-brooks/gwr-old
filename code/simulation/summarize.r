library(plotrix)

vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
params = c('bw', 'sigma2', 'loss.local', 's')

#source('code/matplot.r')

#args = commandArgs(trailingOnly=TRUE)
#cluster = as.integer(args[1])
#cluster = 'NA'
cluster = 56

B = 100
N = 30
coord = seq(0, 1, length.out=N)

#Establish the simulation parameters
tau = rep(c(0.03, 0.1), each=9)
rho = rep(rep(c(0, 0.5, 0.8), each=3), times=2)
sigma.tau = rep(c(0, 0.03, 0.1), times=6)
params = data.frame(tau, rho, sigma.tau)

N = 30
B = list()
settings = 1:18
nsims = 100

coord = seq(0, 1, length.out=N)
B[['(Intercept)']] = rep(0, N**2)
B[['X1']] = as.vector(matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N))
B[['X2']] = rep(0, N**2)
B[['X3']] = rep(0, N**2)
B[['X4']] = rep(0, N**2)
B[['X5']] = rep(0, N**2)

sim.modes = c("gwr", "lars", "enet", "glmnet", "unshrunk.lars", "unshrunk.enet", "unshrunk.glmnet", "oracular")
selection.modes = c("lars", "enet", "glmnet")
file.endings = list(gwr=".gwr.csv", lars=".lars.csv", enet=".enet.csv", glmnet=".glmnet.csv",
                        unshrunk.lars=".unshrunk.lars.csv", unshrunk.enet=".unshrunk.enet.csv",
                        unshrunk.glmnet=".unshrunk.glmnet.csv", oracular=".oracle.csv")

Y.err = list()
X1.err = list()
selection = list()

for (m in sim.modes) {
    Y.err[[m]] = list()
    X1.err[[m]] = list()
}

for (m in selection.modes) {
    selection[[m]] = list()
}

#remove the simulations that failed during their runs.
sims = lapply(1:18, function(x) {as.character(1:100)})
sims[[14]] = sims[[14]][c(-96,-85,-66,-50,-14)]
sims[[11]] = sims[[11]][c(-94,-44)]
sims[[8]] = sims[[8]][-84]


for (setting in settings) {
	cat(paste("Begin setting ", setting, ".\n", sep=""))

    for (m in selection.modes) {
        selection[[m]][[setting]] = list()
    }

    for (m in sim.modes) {
        Y.err[[m]][[setting]] = list()
        X1.err[[m]][[setting]] = list()
    }


    vars = c('(Intercept)', 'X1', 'X2', 'X3', 'X4', 'X5')
    vv = c("X1", 'X2', 'X3', 'X4', 'X5')

    for (k in sims[[setting]]) {
        sim = (setting-1)*100 + as.numeric(k) - 1
        cat(paste("Begin simulation ", sim, ".\n", sep=""))

        #Get the raw data
        filename = paste("output/Data.", cluster, ".", sim, ".csv", sep="")
        raw = read.csv(filename, header=TRUE)

        for (m in sim.modes) {
            filename = paste("output/CoefEstimates.", cluster, ".", sim, file.endings[[m]], sep="")
            estimates = read.csv(filename, header=TRUE)
            colnames(estimates) = vars
            
            if (!is.na(pmatch("unshrunk", m))) {
                fitted = diag(as.matrix(estimates) %*% t(as.matrix(cbind(rep(1,N**2), raw[,c('X1','X2','X3','X4','X5')]))))
            } else {
                filename = paste("output/MiscParams.", cluster, ".", sim, file.endings[[m]], sep="")
                params = read.csv(filename, header=TRUE)
                fitted = params$fitted
            }

            X1.err[[m]][[setting]][[k]] = as.vector(B[['X1']] - estimates$X1)
            Y.err[[m]][[setting]][[k]] = as.vector(raw$Y - fitted)

            if (m %in% selection.modes) {
                for (v in vv) {
                    cat(paste("Begin variable ", v, "\n", sep=""))

                    #Calculate how often each variable was selected for inclusion in the local models
                    col = which(colnames(estimates) == v)
                    if (k==1) {
                        selection[[m]][[setting]][[v]] = as.matrix(ifelse(estimates[,col]==0, 0, 1))
                    } else {
                        selection[[m]][[setting]][[v]] = cbind(selection[[m]][[setting]][[v]], as.matrix(ifelse(estimates[,col]==0, 0, 1)))
                    }
                }
            }
        }
    }
}

