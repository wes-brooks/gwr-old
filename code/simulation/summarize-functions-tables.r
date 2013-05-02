library(xtable)
source("~/git/brooks/code/xtable.printbold.r")
load('~/git/gwr/funcsim2.Rdata')

N = 30
settings = 1:12
nsims = 100

coord = seq(0, 1, length.out=N)

B = list()
for (s in settings) {
    B[[s]] = list()
    ##########################################################
    #Get the true coefficient values for this setting:
    B[[s]][['(Intercept)']] = rep(0, N**2)

    if ((setting-1) %/% 4 == 0) {
        B[[s]][['X1']] = matrix(rep(ifelse(coord<=0.4, 0, ifelse(coord<0.6,5*(coord-0.4),1)), N), N, N)
    } else if ((setting-1) %/% 4 == 1) {
        B[[s]][['X1']] = matrix(rep(coord, N), N, N)
    } else if ((setting-1) %/% 4 == 2) {
        Xmat = matrix(rep(rep(coord, times=N), times=N), N**2, N**2)
        Ymat = matrix(rep(rep(coord, each=N), times=N), N**2, N**2)
        D = (Xmat-t(Xmat))**2 + (Ymat-t(Ymat))**2
        d = D[435,]
        B[[s]][['X1']] = matrix(max(d)-d, N, N)
    }

    B[[s]][['X2']] = rep(0, N**2)
    B[[s]][['X3']] = rep(0, N**2)
    B[[s]][['X4']] = rep(0, N**2)
    B[[s]][['X5']] = rep(0, N**2)
    ##########################################################
}

sim.modes = c("lars", "enet", "glmnet", "unshrunk.lars", "unshrunk.enet", "unshrunk.glmnet", "oracular", "gwr")
selection.modes = c("lars", "enet", "glmnet")
columns = list(lars="LARS", enet="enet", glmnet="glmnet", gwr="gwr",
                        unshrunk.lars=".LARS-U", unshrunk.enet="enet-U",
                        unshrunk.glmnet="glmnet-U", oracular="Oracle")
vars = c('X1', 'X2', 'X3', 'X4', 'X5')

msex = list()
msex.table = list()
msex.fancytable = list()

for (j in 1:length(vars)) {
    msex[[vars[j]]] = list()
    msex.table[[vars[j]]] = list()
    msex.fancytable[[vars[j]]] = list()

    locs = c(30, 228, 435, 643, 871)
    for (l in 1:length(locs)) {
        msex[[vars[j]]][[l]] = list()

        for (m in sim.modes) {
            msex[[vars[j]]][[l]][[m]] = vector()
        }
    }

    for (s in settings) {
        for (l in 1:length(locs)) {
            for (m in sim.modes) {
                msex[[vars[j]]][[l]][[m]] = c(msex[[vars[j]]][[l]][[m]], mean(sapply(X.err[[j]][[m]][[s]], function(x) {x[locs[l]]**2}), na.rm=TRUE))
            }
        }
    }
    
    for (l in 1:length(locs)) {
        msex.table[[vars[j]]][[l]] = as.matrix(sapply(msex[[vars[j]]][[l]], identity))
        msexbold = matrix(FALSE, nrow=dim(msex.table[[vars[j]]][[l]])[1], ncol=dim(msex.table[[vars[j]]][[l]])[2])
        for (i in 1:(dim(msex.table[[vars[j]]][[l]])[1])) {msexbold[i,order(msex.table[[vars[j]]][[l]][i,])[1]] = TRUE}
        msexital = matrix(FALSE, nrow=dim(msex.table[[vars[j]]][[l]])[1], ncol=dim(msex.table[[vars[j]]][[l]])[2])
        for (i in 1:(dim(msex.table[[vars[j]]][[l]])[1])) {msexital[i,order(msex.table[[vars[j]]][[l]][i,])[2]] = TRUE}
        msex.fancytable[[vars[j]]][[l]] = xtable.printbold(xtable(msex.table[[vars[j]]][[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Mean squared error of estimates for $\\beta_" , j, "$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{table:loc", l, "-", vars[j], "-MSEX}", sep="")), which.bold=msexbold, which.ital=msexital, include.rownames=FALSE, hline.after=c(0))
    }
}







msey = list()
msey.table = list()
msey.fancytable = list()

for (l in 1:length(locs)) {
    msey[[l]] = list()

    for (m in sim.modes) {
        msey[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            msey[[l]][[m]] = c(msey[[l]][[m]], mean(sapply(Y.err[[m]][[s]], function(x) {x[locs[l]]**2})))
        }
    }
}

for (l in 1:length(locs)) {
    msey.table[[l]] = as.matrix(sapply(msey[[l]], identity))
    mseybold = matrix(FALSE, nrow=dim(msey.table[[l]])[1], ncol=dim(msey.table[[l]])[2])
    for (i in 1:(dim(msey.table[[l]])[1])) {mseybold[i,order(msey.table[[l]][i,])[1]] = TRUE}
    mseyital = matrix(FALSE, nrow=dim(msey.table[[l]])[1], ncol=dim(msey.table[[l]])[2])
    for (i in 1:(dim(msey.table[[l]])[1])) {mseyital[i,order(msey.table[[l]][i,])[2]] = TRUE}
    msey.fancytable[[l]] = xtable.printbold(xtable(msey.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Mean squared error of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{table:loc", l, "-MSEY}", sep="")), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE, hline.after=c(0))
}







bx = list()
bx.table = list()
bx.fancytable = list()

for (j in 1:length(vars)) {
    bx[[vars[j]]] = list()
    bx.table[[vars[j]]] = list()
    bx.fancytable[[vars[j]]] = list()

    for (l in 1:length(locs)) {
        bx[[vars[j]]][[l]] = list()

        for (m in sim.modes) {
            bx[[vars[j]]][[l]][[m]] = vector()
        }
    }

    for (s in settings) {
        for (l in 1:length(locs)) {
            for (m in sim.modes) {
                bx[[vars[j]]][[l]][[m]] = c(bx[[vars[j]]][[l]][[m]], mean(sapply(X.err[[j]][[m]][[s]], function(x) {-x[locs[l]]})))
            }
        }
    }

    for (l in 1:length(locs)) {
        bx.table[[vars[j]]][[l]] = as.matrix(sapply(bx[[vars[j]]][[l]], identity))
        bxbold = matrix(FALSE, nrow=dim(bx.table[[vars[j]]][[l]])[1], ncol=dim(bx.table[[vars[j]]][[l]])[2])
        for (i in 1:(dim(bx.table[[vars[j]]][[l]])[1])) {bxbold[i,order(abs(bx.table[[vars[j]]][[l]][i,]))[1]] = TRUE}
        bxital = matrix(FALSE, nrow=dim(bx.table[[vars[j]]][[l]])[1], ncol=dim(bx.table[[vars[j]]][[l]])[2])
        for (i in 1:(dim(bx.table[[vars[j]]][[l]])[1])) {bxital[i,order(abs(bx.table[[vars[j]]][[l]][i,]))[2]] = TRUE}
        bx.fancytable[[vars[j]]][[l]] = xtable.printbold(xtable(bx.table[[vars[j]]][[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Bias of estimates for $\\beta_", j, "$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{table:loc", l, "-", vars[[j]], "-BiasX}", sep="")), which.bold=bxbold, which.ital=bxital, include.rownames=FALSE, hline.after=c(0))
    }
}




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
            by[[l]][[m]] = c(by[[l]][[m]], mean(sapply(Y.err[[m]][[s]], function(x) {-x[locs[l]]})))
        }
    }
}

by.table = list()
for (l in 1:length(locs)) {
    by.table[[l]] = as.matrix(sapply(by[[l]], identity))
    bybold = matrix(FALSE, nrow=dim(by.table[[l]])[1], ncol=dim(by.table[[l]])[2])
    for (i in 1:(dim(by.table[[l]])[1])) {bybold[i,order(abs(by.table[[l]][i,]))[1]] = TRUE}
    byital = matrix(FALSE, nrow=dim(by.table[[l]])[1], ncol=dim(by.table[[l]])[2])
    for (i in 1:(dim(by.table[[l]])[1])) {byital[i,order(abs(by.table[[l]][i,]))[2]] = TRUE}
    xtable.printbold(xtable(by.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Bias of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{table:loc", l, "-BiasY}", sep="")), which.bold=bybold, which.ital=byital, include.rownames=FALSE, hline.after=c(0))
}






varx = list()
varx.table = list()
varx.fancytable = list()

for (j in 1:length(vars)) {
    varx[[vars[j]]] = list()
    varx.table[[vars[j]]] = list()
    varx.fancytable[[vars[j]]] = list()

    for (l in 1:length(locs)) {
        varx[[vars[j]]][[l]] = list()

        for (m in sim.modes) {
            varx[[vars[j]]][[l]][[m]] = vector()
        }
    }

    for (s in settings) {
        for (l in 1:length(locs)) {
            for (m in sim.modes) {
                varx[[vars[j]]][[l]][[m]] = c(varx[[vars[j]]][[l]][[m]], var(sapply(X.err[[j]][[m]][[s]], function(x) {x[locs[l]]})))
            }
        }
    }

    for (l in 1:length(locs)) {
        varx.table[[vars[j]]][[l]] = as.matrix(sapply(varx[[vars[j]]][[l]], identity))
        varxbold = matrix(FALSE, nrow=dim(varx.table[[vars[j]]][[l]])[1], ncol=dim(varx.table[[vars[j]]][[l]])[2])
        for (i in 1:(dim(varx.table[[vars[j]]][[l]])[1])) {varxbold[i,order(varx.table[[vars[j]]][[l]][i,])[1]] = TRUE}
        varxital = matrix(FALSE, nrow=dim(varx.table[[vars[j]]][[l]])[1], ncol=dim(varx.table[[vars[j]]][[l]])[2])
        for (i in 1:(dim(varx.table[[vars[j]]][[l]])[1])) {varxital[i,order(varx.table[[vars[j]]][[l]][i,])[2]] = TRUE}
        varx.fancytable[[vars[j]]][[l]] = xtable.printbold(xtable(varx.table[[vars[j]]][[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Variance of estimates for $\\beta_", j, "$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{table:loc", l, "-", vars[[j]], "-varx}", sep="")), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE, hline.after=c(0))
    }
}


vary = list()
vary.table = list()
vary.fancytable = list()

for (l in 1:length(locs)) {
    vary[[l]] = list()

    for (m in sim.modes) {
        vary[[l]][[m]] = vector()
    }
}

for (s in settings) {
    for (l in 1:length(locs)) {
        for (m in sim.modes) {
            vary[[l]][[m]] = c(vary[[l]][[m]], var(sapply(Y.err[[m]][[s]], function(x) {x[locs[l]]})))
        }
    }
}

for (l in 1:length(locs)) {
    vary.table[[l]] = as.matrix(sapply(vary[[l]], identity))
    varybold = matrix(FALSE, nrow=dim(vary.table[[l]])[1], ncol=dim(vary.table[[l]])[2])
    for (i in 1:(dim(vary.table[[l]])[1])) {varybold[i,order(by.table[[l]][i,])[1]] = TRUE}
    varyital = matrix(FALSE, nrow=dim(vary.table[[l]])[1], ncol=dim(vary.table[[l]])[2])
    for (i in 1:(dim(vary.table[[l]])[1])) {varyital[i,order(vary.table[[l]][i,])[2]] = TRUE}
    vary.fancytable[[l]] = xtable.printbold(xtable(vary.table[[l]], digits=3, align=rep('c', length(sim.modes)+1), caption=paste("Variance of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{table:loc", l, "-VarY}", sep="")), which.bold=varybold, which.ital=varyital, include.rownames=FALSE, hline.after=c(0))
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
        comp[[l]][,(m-1)*3+1] = msex[[l]][[sim.modes[m]]] #MSE    
        comp[[l]][,(m-1)*3+2] = bx[[l]][[sim.modes[m]]] #bias
        comp[[l]][,(m-1)*3+3] = varx[[l]][[sim.modes[m]]] #variance    
    }
}


selected = list()
bw = list()

for (v in vars) {
    selected[[v]] = list()

    for (m in selection.modes) {
        selected[[v]][[m]] = list()

        for (l in 1:length(locs)) {
            selected[[v]][[m]][[l]] = vector()
        }        
    }
}

for (m in selection.modes) {
    bw[[m]] = list()

    for (l in 1:length(locs)) {
        bw[[m]][[l]] = vector()
    }        
}

for (s in settings) {
    for (m in selection.modes) {
        for (l in 1:length(locs)) {
            for (v in vars) {
                selected[[v]][[m]][[l]] = c(selected[[v]][[m]][[l]], mean(apply(selection[[m]][[s]][[v]], 2, function(x) {ifelse(x[locs[l]]==0,0,1)})))
            }
            bw[[m]][[l]] = c(bw[[m]][[l]], mean(sapply(bandwidth[[m]][[s]], function(x) {x[locs[l]]})))
        }
    }
}

selection.table = list()
selection.fancytable = list()

for (m in selection.modes) {
    selection.table[[m]] = list()
    selection.fancytable[[m]] = list()

    for (l in 1:length(locs)) {
        selection.table[[m]][[l]] = cbind(bw[[m]][[l]], selected[['X1']][[m]][[l]], rowMeans(sapply(vars[2:5], function(x) {selected[[x]][[m]][[l]]})))
        colnames(selection.table[[m]][[l]]) = c("mean($\\mbox{bw}$)", "$\\beta_1$", "$\\beta_4$ - $\\beta_5$")    
        selection.fancytable[[m]][[l]] = print(xtable(selection.table[[m]][[l]], digits=2, align=c('c','c','c','c'), caption=paste("Selection frequency under ", m, " at location ", l, sep="")), sanitize.colnames.function=function(x){x}, include.rownames=FALSE, hline.after=c(0))
    }
}


comprehensive.selection.table = list()
comprehensive.selection.fancytable = list()

for (l in 1:length(locs)) {
    comprehensive.selection.table[[l]] = list()
    for (m in selection.modes) {
        comprehensive.selection.table[[l]] = c(comprehensive.selection.table[[l]], list(selected[['X1']][[m]][[l]]), list(rowMeans(sapply(vars[2:5], function(x) {selected[[x]][[m]][[l]]}))))
    }
    comprehensive.selection.table[[l]] = as.data.frame(comprehensive.selection.table[[l]])
    colnames(comprehensive.selection.table[[l]]) = rep(c("$\\beta_1$", "$\\beta_4$ - $\\beta_5$"), length(selection.modes))  
    comprehensive.selection.fancytable[[l]] = print(xtable(comprehensive.selection.table[[l]], digits=2, align=c('c', 'c', 'c|', 'c', 'c|', 'c', 'c'), caption=paste("Selection frequency at location ", l, "\\label{table:loc", l, "-selection}", sep="")), sanitize.colnames.function=function(x){x}, include.rownames=FALSE, hline.after=c(0))
}
