library(xtable)
source("~/git/brooks/code/xtable.printbold.r")

N = 30
B = list()
settings = 1:12
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


msex = list()
locs = c(30, 228, 435, 643, 871)
for (l in 1:length(locs)) {
    msex[[l]] = list()

    for (m in sim.modes) {
        msex[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            msex[[l]][[m]] = c(msex[[l]][[m]], mean(sapply(X1.err[[m]][[s]], function(x) {x[locs[l]]**2}), na.rm=TRUE))
        }
    }
}

msex.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (l in 1:length(locs)) {
    msex.table[[l]] = as.matrix(sapply(msex[[l]], identity))
}
msexbold = matrix(FALSE, nrow=nrow(msex.table), ncol=ncol(msex.table))
for (i in 1:nrow(msex.table)) {msexbold[i,order(msex.table[i,])[1]] = TRUE}
msexital = matrix(FALSE, nrow=nrow(msex.table), ncol=ncol(msex.table))
for (i in 1:nrow(msex.table)) {msexital[i,order(msex.table[i,])[2]] = TRUE}
xtable.printbold(xtable(msex.table, digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Mean squared error of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEX}", sep="")), which.bold=msexbold, which.ital=msexital, include.rownames=FALSE, hline.after=c(0))








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

msey.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (l in 1:length(locs)) {
    msey.table = rbind(msey.table, as.matrix(sapply(msey[[l]], identity)))
}
mseybold = matrix(FALSE, nrow=nrow(msey.table), ncol=ncol(msey.table))
for (i in 1:nrow(msey.table)) {mseybold[i,order(msey.table[i,])[1]] = TRUE}
mseyital = matrix(FALSE, nrow=nrow(msey.table), ncol=ncol(msey.table))
for (i in 1:nrow(msey.table)) {mseyital[i,order(msey.table[i,])[2]] = TRUE}
xtable.printbold(xtable(msey.table, digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Mean squared error of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEY}", sep="")), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE, hline.after=c(0))







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

bx.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (l in 1:length(locs)) {
    bx.table = rbind(bx.table, as.matrix(sapply(bx[[l]], identity)))
}

bxbold = matrix(FALSE, nrow=nrow(bx.table), ncol=ncol(bx.table))
for (i in 1:nrow(bx.table)) {bxbold[i,order(bx.table[i,])[1]] = TRUE}
bxital = matrix(FALSE, nrow=nrow(bx.table), ncol=ncol(bx.table))
for (i in 1:nrow(bx.table)) {bxital[i,order(bx.table[i,])[2]] = TRUE}
xtable.printbold(xtable(bx.table, digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Bias of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{BiasX}", sep="")), which.bold=bxbold, which.ital=bxital, include.rownames=FALSE, hline.after=c(0))





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

by.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (l in 1:length(locs)) {
    by.table = rbind(by.table, as.matrix(sapply(by, identity)))
}

bybold = matrix(FALSE, nrow=nrow(by.table), ncol=ncol(by.table))
for (i in 1:nrow(by.table)) {bybold[i,order(by.table[i,])[1]] = TRUE}
byital = matrix(FALSE, nrow=nrow(by.table), ncol=ncol(by.table))
for (i in 1:nrow(by.table)) {byital[i,order(by.table[i,])[2]] = TRUE}
xtable.printbold(xtable(by.table, digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Bias of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{BiasY}", sep="")), which.bold=bybold, which.ital=byital, include.rownames=FALSE, hline.after=c(0))







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

varx.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (l in 1:length(locs)) {
    varx.table = rbind(varx.table, as.matrix(sapply(varx[[l]], identity)))
}

varxbold = matrix(FALSE, nrow=nrow(varx.table), ncol=ncol(varx.table))
for (i in 1:nrow(varx.table)) {varxbold[i,order(varx.table[i,])[1]] = TRUE}
varxital = matrix(FALSE, nrow=nrow(varx.table), ncol=ncol(varx.table))
for (i in 1:nrow(varx.table)) {varxital[i,order(varx.table[i,])[2]] = TRUE}
xtable.printbold(xtable(varx.table, digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Variance of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{VarX}", sep="")), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE, hline.after=c(0))







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

vary.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (l in 1:length(locs)) {
    vary.table = rbind(vary.table, as.matrix(sapply(vary[[l]], identity)))
}

varybold = matrix(FALSE, nrow=nrow(vary.table), ncol=ncol(vary.table))
for (i in 1:nrow(vary.table)) {varybold[i,order(vary.table[i,])[1]] = TRUE}
varyital = matrix(FALSE, nrow=nrow(vary.table), ncol=ncol(vary.table))
for (i in 1:nrow(vary.table)) {varyital[i,order(vary.table[i,])[2]] = TRUE}
xtable.printbold(xtable(vary.table, digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Variance of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{VarY}", sep="")), which.bold=varybold, which.ital=varyital, include.rownames=FALSE, hline.after=c(0))



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

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in selection.modes) {
            for (v in vv) {
                selected[[l]][[m]][[v]] = c(selected[[l]][[m]][[v]], mean(apply(selection[[m]][[s]][[v]], 2, function(x) {ifelse(x[locs[l]]==0,0,1)})))
            }
        }
    }
}

selection.table = matrix(NA, nrow=0, ncol=2*3*selection.modes)
for (j in 1:length(locs)) {
    for (m in selection.modes) {
        selection.table = rbind(selection.table, cbind(selected[[j]][[m]][['X1']], rowMeans(sapply(vv[-1], function(x) {selected[[j]][[m]][[x]]}))))
    }
}
colnames(selection.table) = rep(c("$\\beta_1$", "$\\beta_2$ - $\\beta_5$"), 3*selection.modes)
print(xtable(selection.table, digits=2, align=c('c','c','c'), caption=paste("Selection frequency for the indicated variables at location ", j, sep="")), sanitize.colnames.function=function(x){x}, include.rownames=FALSE, hline.after=c(0))
