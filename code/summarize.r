vars = c('Intercept', 'X1', 'X2', 'X3', 'X4', 'X5', 'Z')
params = c('bw', 'sigma2', 'loss.local', 's', 'sum.weights')

args = commandArgs(trailingOnly=TRUE)
cluster = as.integer(args[1])
#cluster = 25

N = 30
B = list()

coord = seq(0, 1, length.out=N)
B[['Intercept']] = rep(0, N**2)
B[['X1']] = as.vector(matrix(rep(ifelse(coord<=0.5, 0, 2), N), N, N))
B[['X2']] = rep(0, N**2)
B[['X3']] = as.vector(matrix(rep(1-coord, N), N, N))
B[['X4']] = rep(0, N**2)
B[['X5']] = rep(0, N**2)
B[['Z']] = rep(0, N**2)

coverage = list()
selection = list()
indx = 2:7

for (k in indx) {
    #Import our coefficient estimates
    filename = paste("output/CoefEstimates.", cluster, ".", k, ".csv", sep="")
    estimates = read.csv(filename, header=TRUE)
    #estimates = read.table(filename, header=TRUE, sep=' ')

    for (v in vars) {
        #Calculate the coverage of each 95% CI 
        filename = paste("output/", v, ".", cluster, ".", k, ".bootstrap.csv", sep="")
        bootstraps = as.matrix(read.csv(filename, header=FALSE))
        #bootstraps = as.matrix(read.table(filename, header=TRUE, sep=' '))
        CI = t(apply(bootstraps, 1, function(x) {sort(x)[c(4, 98)]}))
        
        if (k==1) {
            coverage[[v]] = as.matrix(ifelse(B[[v]] < CI[,1] | B[[v]] > CI[,2],0,1))
        } else {
            coverage[[v]] = cbind(coverage[[v]], as.matrix(ifelse(B[[v]] < CI[,1] | B[[v]] > CI[,2],0,1)))
        }


        #Calculate how often each variable was selected for inclusion in the local models
        col = which(colnames(estimates) == v)
        if (k==1) {
            selection[[v]] = as.matrix(ifelse(estimates[,col]==0, 0, 1))
        } else {
            selection[[v]] = cbind(selection[[v]], as.matrix(ifelse(estimates[,col]==0, 0, 1)))
        }        
    }
}


for (v in vars) {
    cc = apply(coverage[[v]], 1, sum) / length(indx)
    ss = apply(selection[[v]], 1, sum) / length(indx)

    
}
