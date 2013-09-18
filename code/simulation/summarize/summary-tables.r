library(xtable)
source("~/git/brooks/code/xtable.printbold.r")

#Load the simulation results
args = commandArgs(trailingOnly=TRUE)
root = paste(as.character(args[1]), "/", sep="")
cluster = as.integer(args[2])
outdir = as.character(args[3])
load(paste(root, "sim-", cluster, ".Rdata", sep=""))

N = 30
B = list()
settings = 1:12
nsims = 100

#coord = seq(0, 1, length.out=N)
#B[['(Intercept)']] = rep(0, N**2)
#B[['X1']] = as.vector(matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N))
#B[['X2']] = rep(0, N**2)
#B[['X3']] = rep(0, N**2)
#B[['X4']] = rep(0, N**2)
#B[['X5']] = rep(0, N**2)

#functions = c('step', 'gradient', 'parabola')
#sim.modes = c("lars", "enet", "glmnet", "unshrunk.lars", "unshrunk.enet", "unshrunk.glmnet", "oracular", "gwr")
#sim.modes.output = c("lars", "enet", "glmnet", "u.lars", "u.enet", "u.glmnet", "oracular", "gwr")
#selection.modes = c("lars", "enet", "glmnet")
#columns = list(lars="LARS", enet="enet", glmnet="glmnet", gwr="gwr",
#                        unshrunk.lars=".LARS-U", unshrunk.enet="enet-U",
#                        unshrunk.glmnet="glmnet-U", oracular="Oracle")
#groupings = list('1'=c(1,2,3,4), '2'=c(5,6,7,8), '3'=c(9,10,11,12))


#no lars:
functions = c('step', 'gradient', 'parabola')
sim.modes = c("enet", "glmnet", "unshrunk.enet", "unshrunk.glmnet", "oracular", "gwr")
sim.modes.output = c("enet", "glmnet", "u.enet", "u.glmnet", "oracular", "gwr")
selection.modes = c("enet", "glmnet")
columns = list(enet="enet", glmnet="glmnet", gwr="gwr",
                        unshrunk.enet="enet-U",
                        unshrunk.glmnet="glmnet-U", oracular="Oracle")
groupings = list('1'=c(1,2,3,4), '2'=c(5,6,7,8), '3'=c(9,10,11,12))


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
            msex[[l]][[m]] = c(msex[[l]][[m]], mean(sapply(X.err[[1]][[m]][[s]], function(x) {x[locs[l]]**2}), na.rm=TRUE))
        }
    }
}

msex.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (g in groupings) {
    for (l in 1:length(locs)) {
        msex.table = rbind(msex.table, as.matrix(sapply(sim.modes, function(x) {msex[[l]][[x]][g]})))
    }
}
msexbold = matrix(FALSE, nrow=nrow(msex.table), ncol=ncol(msex.table)+2)
msexital = matrix(FALSE, nrow=nrow(msex.table), ncol=ncol(msex.table)+2)
for (i in 1:nrow(msex.table)) {
    row = msex.table[i,]
    best = sort(unique(row))
    msexbold[i,-(1:2)][row==best[1]] = TRUE
    if (sum(row==best[1], na.rm=TRUE)==1) {msexital[i,-(1:2)][row==best[2]] = TRUE}
}


msex.table = round(msex.table, 3)
msex.table = cbind(rep(NA,nrow(msex.table)), rep(NA,nrow(msex.table)), msex.table)
colnames(msex.table) = c('function','location',sim.modes.output)
msex = xtable.printbold(xtable(msex.table, digits=3, align=rep('c', length(sim.modes)+3), caption="Mean squared error of $\\hat{\\beta_1}$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX}"), which.bold=msexbold, which.ital=msexital, include.rownames=FALSE, hline.after=c(0))

sink(paste(outdir, "/msex.tex", sep=""))
print(msex)
sink()






msey = list()
for (l in 1:length(locs)) {
    msey[[l]] = list()

    for (m in sim.modes) {
        msey[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            msey[[l]][[m]] = c(msey[[l]][[m]], mean(sapply(Y.err[[m]][[s]], function(x) {x[locs[l]]**2}), na.rm=TRUE))
        }
    }
}

msey.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (g in groupings) {
    for (l in 1:length(locs)) {
        msey.table = rbind(msey.table, as.matrix(sapply(sim.modes, function(x) {msey[[l]][[x]][g]})))
    }
}

mseybold = matrix(FALSE, nrow=nrow(msey.table), ncol=ncol(msey.table)+2)
mseyital = matrix(FALSE, nrow=nrow(msey.table), ncol=ncol(msey.table)+2)
for (i in 1:nrow(msey.table)) {
    row = msey.table[i,]
    best = sort(unique(row))
    mseybold[i,-(1:2)][row==best[1]] = TRUE
    if (sum(row==best[1], na.rm=TRUE)==1) {mseyital[i,-(1:2)][row==best[2]] = TRUE}
}

msey.table = round(msey.table, 3)
msey.table = cbind(rep(NA,nrow(msey.table)), rep(NA,nrow(msey.table)), msey.table)
colnames(msey.table) = c('function','location',sim.modes.output)
xtable.printbold(xtable(msey.table, digits=3, align=rep('c', length(sim.modes)+3), caption="Mean squared error of $\\hat{Y}$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEY}"), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE, hline.after=c(0))







bx = list()
for (l in 1:length(locs)) {
    bx[[l]] = list()

    for (m in sim.modes) {
        bx[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            bx[[l]][[m]] = c(bx[[l]][[m]], mean(sapply(X.err[[1]][[m]][[s]], function(x) {-x[locs[l]]}), na.rm=TRUE))
        }
    }
}

bx.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (g in groupings) {
    for (l in 1:length(locs)) {
        bx.table = rbind(bx.table, as.matrix(sapply(sim.modes, function(x) {bx[[l]][[x]][g]})))
    }
}

bxbold = matrix(FALSE, nrow=nrow(bx.table), ncol=ncol(bx.table)+2)
bxital = matrix(FALSE, nrow=nrow(bx.table), ncol=ncol(bx.table)+2)
for (i in 1:nrow(bx.table)) {
    row = bx.table[i,]
    best = sort(unique(abs(row)))
    bxbold[i,-(1:2)][abs(row)==best[1]] = TRUE
    if (sum(abs(row)==best[1], na.rm=TRUE)==1) {bxital[i,-(1:2)][abs(row)==best[2]] = TRUE}
}

bx.table = round(bx.table, 3)
bx.table = cbind(rep(NA,nrow(bx.table)), rep(NA,nrow(bx.table)),bx.table)
colnames(bx.table) = c('function','location',sim.modes.output)
xtable.printbold(xtable(bx.table, digits=3, align=rep('c', length(sim.modes)+3), caption="Bias of $\\hat{\\beta_1}$ (\\textbf{minimum}, \\emph{next best}).\\label{BiasX}"), which.bold=bxbold, which.ital=bxital, include.rownames=FALSE, hline.after=c(0))





by = list()
for (l in 1:length(locs)) {
    by[[l]] = list()

    for (m in sim.modes) {
        by[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            by[[l]][[m]] = c(by[[l]][[m]], mean(sapply(Y.err[[m]][[s]], function(x) {-x[locs[l]]}), na.rm=TRUE))
        }
    }
}

by.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (g in groupings) {
    for (l in 1:length(locs)) {
        by.table = rbind(by.table, as.matrix(sapply(sim.modes, function(x) {by[[l]][[x]][g]})))
    }
}

bybold = matrix(FALSE, nrow=nrow(by.table), ncol=ncol(by.table)+2)
byital = matrix(FALSE, nrow=nrow(by.table), ncol=ncol(by.table)+2)
for (i in 1:nrow(by.table)) {
    row = by.table[i,]
    best = sort(unique(abs(row)))
    bybold[i,-(1:2)][abs(row)==best[1]] = TRUE
    if (sum(abs(row)==best[1], na.rm=TRUE)==1) {byital[i,-(1:2)][abs(row)==best[2]] = TRUE}
}

by.table = round(by.table, 3)
by.table = cbind(rep(NA,nrow(by.table)), rep(NA,nrow(by.table)),by.table)
colnames(by.table) = c('function','location',sim.modes.output)
xtable.printbold(xtable(by.table, digits=3, align=rep('c', length(sim.modes)+3), caption="Bias of $\\hat{Y}$ (\\textbf{minimum}, \\emph{next best}).\\label{BiasY}"), which.bold=bybold, which.ital=byital, include.rownames=FALSE, hline.after=c(0))







varx = list()
for (l in 1:length(locs)) {
    varx[[l]] = list()

    for (m in sim.modes) {
        varx[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            varx[[l]][[m]] = c(varx[[l]][[m]], var(sapply(X.err[[1]][[m]][[s]], function(x) {x[locs[l]]}), na.rm=TRUE))
        }
    }
}

varx.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (g in groupings) {
    for (l in 1:length(locs)) {
        varx.table = rbind(varx.table, as.matrix(sapply(sim.modes, function(x) {varx[[l]][[x]][g]})))
    }
}

varxbold = matrix(FALSE, nrow=nrow(varx.table), ncol=ncol(varx.table)+2)
varxital = matrix(FALSE, nrow=nrow(varx.table), ncol=ncol(varx.table)+2)
for (i in 1:nrow(varx.table)) {
    row = varx.table[i,]
    best = sort(unique(row))
    varxbold[i,-(1:2)][row==best[1]] = TRUE
    if (sum(row==best[1], na.rm=TRUE)==1) {varxital[i,-(1:2)][row==best[2]] = TRUE}
}

varx.table = round(varx.table, 3)
varx.table = cbind(rep(NA,nrow(varx.table)), rep(NA,nrow(varx.table)), varx.table)
colnames(varx.table) = c('function','location',sim.modes.output)
xtable.printbold(xtable(varx.table, digits=3, align=rep('c', length(sim.modes)+3), caption="Variance of $\\hat{\\beta_1}$ (\\textbf{minimum}, \\emph{next best}).\\label{VarX}"), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE, hline.after=c(0))







vary = list()
for (l in 1:length(locs)) {
    vary[[l]] = list()

    for (m in sim.modes) {
        vary[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            vary[[l]][[m]] = c(vary[[l]][[m]], var(sapply(Y.err[[m]][[s]], function(x) {x[locs[l]]}), na.rm=TRUE))
        }
    }
}

vary.table = matrix(NA, nrow=0, ncol=length(sim.modes))
for (g in groupings) {
    for (l in 1:length(locs)) {
        vary.table = rbind(vary.table, as.matrix(sapply(sim.modes, function(x) {vary[[l]][[x]][g]})))
    }
}

varybold = matrix(FALSE, nrow=nrow(vary.table), ncol=ncol(vary.table)+2)
varyital = matrix(FALSE, nrow=nrow(vary.table), ncol=ncol(vary.table)+2)
for (i in 1:nrow(vary.table)) {
    row = vary.table[i,]
    best = sort(unique(row))
    varybold[i,-(1:2)][row==best[1]] = TRUE
    if (sum(row==best[1], na.rm=TRUE)==1) {varyital[i,-(1:2)][row==best[2]] = TRUE}
}

vary.table = round(vary.table, 3)
vary.table = cbind(rep(NA,nrow(vary.table)), rep(NA,nrow(vary.table)), vary.table)
colnames(vary.table) = c('function','location',sim.modes.output)
xtable.printbold(xtable(vary.table, digits=3, align=rep('c', length(sim.modes)+3), caption="Variance of $\\hat{Y}$ (\\textbf{minimum}, \\emph{next best}).\\label{VarY}"), which.bold=varybold, which.ital=varyital, include.rownames=FALSE, hline.after=c(0))



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

selection.table = matrix(NA, nrow=length(groupings[[1]])*length(locs), ncol=0)
for (j in 1:3) {    
    for (m in selection.modes) {
        selection.block = matrix(NA, nrow=0, ncol=2)
        for (l in 1:length(locs)) {
            k = length(groupings[[j]])
            for (i in 1:k) {
                selection.block = rbind(selection.block, c(selected[[l]][[m]][['X1']][(j-1)*k + i], mean(sapply(vv[-1], function(x) {selected[[l]][[m]][[x]][(j-1)*k + i]}))))
            }
        }
        selection.table = cbind(selection.table, selection.block)
    }
}

selection.table = round(selection.table, 2)
selection.table = cbind(rep(NA,nrow(selection.table)), selection.table)
colnames(selection.table) = c('location', rep(c("$\\beta_1$", "$\\beta_2$ - $\\beta_5$"), length(functions)*length(selection.modes)))
print(xtable(selection.table, digits=2, align=rep('c',ncol(selection.table)+1), caption="Selection frequency for the indicated variables"), sanitize.colnames.function=function(x){x}, include.rownames=FALSE, hline.after=c(0))
