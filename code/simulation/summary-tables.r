library(xtable)
source("~/git/brooks/code/xtable.printbold.r")

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

sim.modes = c("lars", "enet", "glmnet", "unshrunk.lars", "unshrunk.enet", "unshrunk.glmnet", "oracular", "gwr")
selection.modes = c("lars", "enet", "glmnet")
columns = list(lars="LARS", enet="enet", glmnet="glmnet", gwr="gwr",
                        unshrunk.lars=".LARS-U", unshrunk.enet="enet-U",
                        unshrunk.glmnet="glmnet-U", oracular="Oracle")


mse = list()
locs = c(30, 228, 435, 643, 871)
for (l in 1:length(locs)) {
    mse[[l]] = list()

    for (m in sim.modes) {
        mse[[l]][[m]] = vector()
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            mse[[l]][[m]] = c(mse[[l]][[m]], mean(sapply(X1.err[[m]][[s]], function(x) {x[locs[l]]**2}), na.rm=TRUE))
        }
    }
}

mse.table = list()
for (l in 1:length(locs)) {
    mse.table[[l]] = as.matrix(sapply(mse[[l]], identity))
    msebold = matrix(FALSE, nrow=dim(mse.table[[l]])[1], ncol=dim(mse.table[[l]])[2])
    for (i in 1:(dim(mse.table[[l]])[1])) {msebold[i,order(mse.table[[l]][i,])[1]] = TRUE}
    mseital = matrix(FALSE, nrow=dim(mse.table[[l]])[1], ncol=dim(mse.table[[l]])[2])
    for (i in 1:(dim(mse.table[[l]])[1])) {mseital[i,order(mse.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(mse.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Mean squared error of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEX}", sep="")), which.bold=msebold, which.ital=mseital, include.rownames=FALSE, hline.after=c(0))
}








msey = list()
for (l in 1:length(locs)) {
    msey[[l]] = list()

    for (m in sim.modes) {
        msey[[l]][[m]] = vector()
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            msey[[l]][[m]] = c(msey[[l]][[m]], mean(sapply(Y.err[[m]][[s]], function(x) {x[locs[l]]**2})))
        }
    }
}

msey.table = list()
for (l in 1:length(locs)) {
    msey.table[[l]] = as.matrix(sapply(msey[[l]], identity))
    mseybold = matrix(FALSE, nrow=dim(msey.table[[l]])[1], ncol=dim(msey.table[[l]])[2])
    for (i in 1:(dim(msey.table[[l]])[1])) {mseybold[i,order(msey.table[[l]][i,])[1]] = TRUE}
    mseyital = matrix(FALSE, nrow=dim(msey.table[[l]])[1], ncol=dim(msey.table[[l]])[2])
    for (i in 1:(dim(msey.table[[l]])[1])) {mseyital[i,order(msey.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(msey.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Mean squared error of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEY}", sep="")), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE, hline.after=c(0))
}







bx = list()
for (l in 1:length(locs)) {
    bx[[l]] = list()

    for (m in sim.modes) {
        bx[[l]][[m]] = vector()
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            bx[[l]][[m]] = c(bx[[l]][[m]], mean(sapply(X1.err[[m]][[s]], function(x) {-x[locs[l]]})))
        }
    }
}

bx.table = list()
for (l in 1:length(locs)) {
    bx.table[[l]] = as.matrix(sapply(bx[[l]], identity))
    bxbold = matrix(FALSE, nrow=dim(bx.table[[l]])[1], ncol=dim(bx.table[[l]])[2])
    for (i in 1:(dim(bx.table[[l]])[1])) {bxbold[i,order(bx.table[[l]][i,])[1]] = TRUE}
    bxital = matrix(FALSE, nrow=dim(bx.table[[l]])[1], ncol=dim(bx.table[[l]])[2])
    for (i in 1:(dim(bx.table[[l]])[1])) {bxital[i,order(bx.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(bx.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Bias of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{BiasX}", sep="")), which.bold=bxbold, which.ital=bxital, include.rownames=FALSE, hline.after=c(0))
}





by = list()
for (l in 1:length(locs)) {
    by[[l]] = list()

    for (m in sim.modes) {
        by[[l]][[m]] = vector()
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            by[[l]][[m]] = c(by[[l]][[m]], mean(sapply(Y.err[[m]][[s]], function(x) {-x[locs[l]]})))
        }
    }
}

by.table = list()
for (l in 1:length(locs)) {
    by.table[[l]] = as.matrix(sapply(by[[l]], identity))
    bybold = matrix(FALSE, nrow=dim(by.table[[l]])[1], ncol=dim(by.table[[l]])[2])
    for (i in 1:(dim(by.table[[l]])[1])) {bybold[i,order(by.table[[l]][i,])[1]] = TRUE}
    byital = matrix(FALSE, nrow=dim(by.table[[l]])[1], ncol=dim(by.table[[l]])[2])
    for (i in 1:(dim(by.table[[l]])[1])) {byital[i,order(by.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(by.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Bias of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{BiasY}", sep="")), which.bold=bybold, which.ital=byital, include.rownames=FALSE, hline.after=c(0))
}







varx = list()
for (l in 1:length(locs)) {
    varx[[l]] = list()

    for (m in sim.modes) {
        varx[[l]][[m]] = vector()
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            varx[[l]][[m]] = c(varx[[l]][[m]], var(sapply(X1.err[[m]][[s]], function(x) {x[locs[l]]})))
        }
    }
}

varx.table = list()
for (l in 1:length(locs)) {
    varx.table[[l]] = as.matrix(sapply(varx[[l]], identity))
    varxbold = matrix(FALSE, nrow=dim(varx.table[[l]])[1], ncol=dim(varx.table[[l]])[2])
    for (i in 1:(dim(varx.table[[l]])[1])) {varxbold[i,order(by.table[[l]][i,])[1]] = TRUE}
    varxital = matrix(FALSE, nrow=dim(varx.table[[l]])[1], ncol=dim(varx.table[[l]])[2])
    for (i in 1:(dim(varx.table[[l]])[1])) {varxital[i,order(varx.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(varx.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Variance of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{VarX}", sep="")), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE, hline.after=c(0))
}







vary = list()
for (l in 1:length(locs)) {
    vary[[l]] = list()

    for (m in sim.modes) {
        vary[[l]][[m]] = vector()
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            vary[[l]][[m]] = c(vary[[l]][[m]], var(sapply(X1.err[[m]][[s]], function(x) {x[locs[l]]})))
        }
    }
}

vary.table = list()
for (l in 1:length(locs)) {
    vary.table[[l]] = as.matrix(sapply(vary[[l]], identity))
    varybold = matrix(FALSE, nrow=dim(vary.table[[l]])[1], ncol=dim(vary.table[[l]])[2])
    for (i in 1:(dim(vary.table[[l]])[1])) {varybold[i,order(by.table[[l]][i,])[1]] = TRUE}
    varyital = matrix(FALSE, nrow=dim(vary.table[[l]])[1], ncol=dim(vary.table[[l]])[2])
    for (i in 1:(dim(vary.table[[l]])[1])) {varyital[i,order(vary.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(vary.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Variance of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{VarY}", sep="")), which.bold=varybold, which.ital=varyital, include.rownames=FALSE, hline.after=c(0))
}





comprehensive.table = list()
for (l in 1:length(locs)) {
    for (m in sim.modes) {
        comprehensive.table[[l]] = cbind(mse.loc[[l]][['GWL']], bx.loc[[l]][['GWL']], varx.loc[[l]][['GWL']], mse.loc[[l]][['unshrunk']], bx.loc[[l]][['unshrunk']], varx.loc[[l]][['unshrunk']], mse.loc[[l]][['precon']], bx.loc[[l]][['precon']], varx.loc[[l]][['precon']], mse.loc[[l]][['unshrunk-precon']], bx.loc[[l]][['unshrunk-precon']], varx.loc[[l]][['unshrunk-precon']], mse.loc[[l]][['oracular']], bx.loc[[l]][['oracular']], varx.loc[[l]][['oracular']])
        colnames(comprehensive.table[[l]]) = c(rep("GWL", 3), rep("GWL-U", 3), rep("GWL-P", 3), rep("GWL-P-U", 3), rep("Oracle", 3))
        compbold = matrix(FALSE, nrow=dim(comprehensive.loc.table2[[l]])[1], ncol=dim(comprehensive.loc.table2[[l]])[2])
        for (j in 1:3) {
            cols = dim(comprehensive.loc.table2[[l]])[2]
            indx = which((1:cols - 1) %% 3 + 1 == j)        
            for (i in 1:(dim(comprehensive.loc.table2[[l]])[1])) {
                minloc = indx[order(abs(comprehensive.loc.table2[[l]][i,indx]))[1]]
                compbold[i,minloc] = TRUE
            }
        }
    }
    compital = matrix(FALSE, nrow=dim(comprehensive.loc.table2[[l]])[1], ncol=dim(comprehensive.loc.table2[[l]])[2])
    for (j in 1:3) {
        cols = dim(comprehensive.loc.table2[[l]])[2]
        indx = which((1:cols - 1) %% 3 + 1 == j)        
        for (i in 1:(dim(comprehensive.loc.table2[[l]])[1])) {
            minloc = indx[order(abs(comprehensive.loc.table2[[l]][i,indx]))[2]]
            compital[i,minloc] = TRUE
        }
    }
    xtable.printbold(xtable(comprehensive.loc.table2[[l]], digits=3, align=c('c','c','c','c','c','c','c','c','c','c','c','c','c','c','c','c'), caption=paste("MSE, bias, and variance of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).", sep="")), include.rownames=FALSE, hline.after=c(0), which.bold=compbold, which.ital=compital)
    #print(xtable(comprehensive.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c','c','c','c','c'), caption=paste("MSE, bias, and variance of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).", sep="")), include.rownames=FALSE, hline.after=c(0))
}

comp = list()
for (l in 1:length(locs)) {
    comp[[l]] = matrix(NA, nrow=length(settings), ncol=length(sim.modes)*3)
    for (m in 1:length(sim.modes)) {
        comp[[l]][,(m-1)*3+1] = mse[[l]][[sim.modes[m]]] #MSE    
        comp[[l]][,(m-1)*3+2] = bx[[l]][[sim.modes[m]]] #bias
        comp[[l]][,(m-1)*3+3] = varx[[l]][[sim.modes[m]]] #variance    
    }
}


vv = c('X1', 'X2', 'X3', 'X4', 'X5')
selected = list()

for (l in 1:length(locs)) {
    selected[[l]] = list()

    for (m in selection.modes) {
        selected[[l]][[m]] = list()
        for (v in vv) {
            selected[[l]][[m]][[v]] = vector()
        }        
    }
}

for (s in 1:18) {
    for (l in 1:length(locs)) {
        for (m in selection.modes) {
            for (v in vv) {
                selected[[l]][[m]][[v]] = c(selected[[l]][[m]][[v]], mean(apply(selection[[m]][[s]][[v]], 2, function(x) {ifelse(x[locs[l]]==0,0,1)})))
            }
        }
    }
}

selection.table = list()
for (j in 1:length(locs)) {
    for (m in selection.modes) {
        selection.table[[j]] = cbind(selected[[j]][[m]][['X1']], rowMeans(sapply(vv[-1], function(x) {selected[[j]][[m]][[x]]})))
        colnames(selection.table[[j]]) = c("$\\beta_1$", "$\\beta_2$ - $\\beta_5$")    
        print(xtable(selection.table[[j]], digits=2, align=c('c','c','c'), caption=paste("Selection frequency for the indicated variables at location ", j, sep="")), sanitize.colnames.function=function(x){x}, include.rownames=FALSE, hline.after=c(0))
    }
}
