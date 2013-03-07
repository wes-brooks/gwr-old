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


mse = vector()
mse.oracular = vector()
mse.precon = vector()
mse.unshrunk = vector()
mse.unshrunk.precon = vector()

locs = c(30, 318, 450, 613, 871)
mse.loc = list()
for (l in 1:length(locs)) {    
    mse.loc[[l]] = list()
    mse.loc[[l]][['GWL']] = vector()
    mse.loc[[l]][['oracular']] = vector()
    mse.loc[[l]][['precon']] = vector()
    mse.loc[[l]][['unshrunk']] = vector()
    mse.loc[[l]][['unshrunk-precon']] = vector()
}

mse.nonzero = vector()
mse.oracular.nonzero = vector()
mse.precon.nonzero = vector()
mse.unshrunk.nonzero = vector()
mse.unshrunk.precon.nonzero = vector()

for (s in 1:18) {
	mse = c(mse, mean(sapply(X1.err[[s]], function(x) {mean(x**2)})))
	mse.oracular = c(mse.oracular, mean(sapply(X1.err.oracular[[s]], function(x) {mean(x**2)})))
	mse.precon = c(mse.precon, mean(sapply(X1.err.precon[[s]], function(x) {mean(x**2)})))
	mse.unshrunk = c(mse.unshrunk, mean(sapply(X1.err.unshrunk[[s]], function(x) {mean(x**2)})))
	mse.unshrunk.precon = c(mse.unshrunk.precon, mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {mean(x**2)})))

    for (l in 1:length(locs)) {
        mse.loc[[l]][['GWL']] = c(mse.loc[[l]][['GWL']], mean(sapply(X1.err[[s]], function(x) {x[locs[l]]**2})))
        mse.loc[[l]][['oracular']] = c(mse.loc[[l]][['oracular']], mean(sapply(X1.err.oracular[[s]], function(x) {x[locs[l]]**2})))
        mse.loc[[l]][['precon']] = c(mse.loc[[l]][['precon']], mean(sapply(X1.err.precon[[s]], function(x) {x[locs[l]]**2})))
        mse.loc[[l]][['unshrunk']] = c(mse.loc[[l]][['unshrunk']], mean(sapply(X1.err.unshrunk[[s]], function(x) {x[locs[l]]**2})))
        mse.loc[[l]][['unshrunk-precon']] = c(mse.loc[[l]][['unshrunk-precon']], mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[locs[l]]**2})))
    }

	mse.nonzero = c(mse.nonzero, mean(sapply(X1.err[[s]], function(x) {mean(x[B[['X1']]!=0]**2)}), na.rm=TRUE))
    mse.oracular.nonzero = c(mse.oracular.nonzero, mean(sapply(X1.err.oracular[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
	mse.precon.nonzero = c(mse.precon.nonzero, mean(sapply(X1.err.precon[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
	mse.unshrunk.nonzero = c(mse.unshrunk.nonzero, mean(sapply(X1.err.unshrunk[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
	mse.unshrunk.precon.nonzero = c(mse.unshrunk.precon.nonzero, mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
}
mse2.table = cbind(mse, mse.unshrunk, mse.oracular)
colnames(mse2.table) = c("GWL", "GWL-U", "Oracle")
mse2bold = matrix(FALSE, nrow=dim(mse2.table)[1], ncol=dim(mse2.table)[2])
for (i in 1:(dim(mse2.table)[1])) {mse2bold[i,order(mse2.table[i,])[1]] = TRUE}
mse2ital = matrix(FALSE, nrow=dim(mse2.table)[1], ncol=dim(mse2.table)[2])
for (i in 1:(dim(mse2.table)[1])) {mse2ital[i,order(mse2.table[i,])[2]] = TRUE}
xtable.printbold(xtable(mse2.table, digits=3, align=c('c','c','c','c'), caption="Mean squared error of estimates for $\\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX}"), which.bold=mse2bold, which.ital=mse2ital, include.rownames=FALSE, hline.after=c(0))

mse2.nonzero.table = cbind(mse.nonzero, mse.unshrunk.nonzero, mse.oracular.nonzero)
colnames(mse2.nonzero.table) = c("GWL", "GWL-U", "Oracle")
mse2bold.nonzero = matrix(FALSE, nrow=dim(mse2.nonzero.table)[1], ncol=dim(mse2.nonzero.table)[2])
for (i in 1:(dim(mse2.nonzero.table)[1])) {mse2bold.nonzero[i,order(mse2.nonzero.table[i,])[1]] = TRUE}
mse2ital.nonzero = matrix(FALSE, nrow=dim(mse2.nonzero.table)[1], ncol=dim(mse2.nonzero.table)[2])
for (i in 1:(dim(mse2.nonzero.table)[1])) {mse2ital.nonzero[i,order(mse2.nonzero.table[i,])[2]] = TRUE}
xtable.printbold(xtable(mse2.nonzero.table, digits=3, align=c('c','c','c','c'), caption="Mean squared error of estimates for $\\beta_1$ at locations where $\\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}"), which.bold=mse2bold.nonzero, which.ital=mse2ital.nonzero, include.rownames=FALSE, hline.after=c(0))


mse.table = cbind(mse, mse.unshrunk, mse.precon, mse.unshrunk.precon, mse.oracular)
colnames(mse.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
msebold = matrix(FALSE, nrow=dim(mse.table)[1], ncol=dim(mse.table)[2])
for (i in 1:(dim(mse.table)[1])) {msebold[i,order(mse.table[i,])[1]] = TRUE}
mseital = matrix(FALSE, nrow=dim(mse.table)[1], ncol=dim(mse.table)[2])
for (i in 1:(dim(mse.table)[1])) {mseital[i,order(mse.table[i,])[2]] = TRUE}
xtable.printbold(xtable(mse.table, digits=3, align=c('c','c','c','c','c','c'), caption="Mean squared error of estimates for $\\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}"), which.bold=msebold, which.ital=mseital, include.rownames=FALSE, hline.after=c(0))

mse.nonzero.table = cbind(mse.nonzero, mse.unshrunk.nonzero, mse.precon.nonzero, mse.unshrunk.precon.nonzero, mse.oracular.nonzero)
colnames(mse.nonzero.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
msebold.nonzero = matrix(FALSE, nrow=dim(mse.nonzero.table)[1], ncol=dim(mse.nonzero.table)[2])
for (i in 1:(dim(mse.nonzero.table)[1])) {msebold.nonzero[i,order(mse.nonzero.table[i,])[1]] = TRUE}
mseital.nonzero = matrix(FALSE, nrow=dim(mse.nonzero.table)[1], ncol=dim(mse.nonzero.table)[2])
for (i in 1:(dim(mse.nonzero.table)[1])) {mseital.nonzero[i,order(mse.nonzero.table[i,])[2]] = TRUE}
xtable.printbold(xtable(mse.nonzero.table, digits=3, align=c('c','c','c','c','c','c'), caption="Mean squared error of estimates for $\\beta_1$ at locations where $\\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}"), which.bold=msebold.nonzero, which.ital=mseital.nonzero, include.rownames=FALSE, hline.after=c(0))


mse.loc.table = list()
for (l in 1:length(locs)) {
    mse.loc.table[[l]] = cbind(mse.loc[[l]][['GWL']], mse.loc[[l]][['unshrunk']], mse.loc[[l]][['precon']], mse.loc[[l]][['unshrunk-precon']], mse.loc[[l]][['oracular']])
    colnames(mse.loc.table[[l]]) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
    msebold = matrix(FALSE, nrow=dim(mse.loc.table[[l]])[1], ncol=dim(mse.loc.table[[l]])[2])
    for (i in 1:(dim(mse.loc.table[[l]])[1])) {msebold[i,order(mse.loc.table[[l]][i,])[1]] = TRUE}
    mseital = matrix(FALSE, nrow=dim(mse.loc.table[[l]])[1], ncol=dim(mse.loc.table[[l]])[2])
    for (i in 1:(dim(mse.loc.table[[l]])[1])) {mseital[i,order(mse.loc.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(mse.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c'), caption=paste("Mean squared error of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}", sep="")), which.bold=msebold, which.ital=mseital, include.rownames=FALSE, hline.after=c(0))
}




locs = c(30, 318, 450, 613, 871)
msey.loc = list()
for (l in 1:length(locs)) {    
    msey.loc[[l]] = list()
    msey.loc[[l]][['GWL']] = vector()
    msey.loc[[l]][['oracular']] = vector()
    msey.loc[[l]][['precon']] = vector()
    msey.loc[[l]][['unshrunk']] = vector()
    msey.loc[[l]][['unshrunk-precon']] = vector()
}

msey = vector()
msey.oracular = vector()
msey.precon = vector()
msey.unshrunk = vector()
msey.unshrunk.precon = vector()
for (s in 1:18) {
	msey = c(msey, mean(sapply(Y.err[[s]], function(x) {mean(x**2)})))
	msey.oracular = c(msey.oracular, mean(sapply(Y.err.oracular[[s]], function(x) {mean(x**2)})))
	msey.precon = c(msey.precon, mean(sapply(Y.err.precon[[s]], function(x) {mean(x**2)})))
	msey.unshrunk = c(msey.unshrunk, mean(sapply(Y.err.unshrunk[[s]], function(x) {mean(x**2)})))
	msey.unshrunk.precon = c(msey.unshrunk.precon, mean(sapply(Y.err.unshrunk.precon[[s]], function(x) {mean(x**2)})))

    for (l in 1:length(locs)) {
        msey.loc[[l]][['GWL']] = c(msey.loc[[l]][['GWL']], mean(sapply(Y.err[[s]], function(x) {x[locs[l]]**2})))
        msey.loc[[l]][['oracular']] = c(msey.loc[[l]][['oracular']], mean(sapply(Y.err.oracular[[s]], function(x) {x[locs[l]]**2})))
        msey.loc[[l]][['precon']] = c(msey.loc[[l]][['precon']], mean(sapply(Y.err.precon[[s]], function(x) {x[locs[l]]**2})))
        msey.loc[[l]][['unshrunk']] = c(msey.loc[[l]][['unshrunk']], mean(sapply(Y.err.unshrunk[[s]], function(x) {x[locs[l]]**2})))
        msey.loc[[l]][['unshrunk-precon']] = c(msey.loc[[l]][['unshrunk-precon']], mean(sapply(Y.err.unshrunk.precon[[s]], function(x) {x[locs[l]]**2})))
    }
}
msey2.table = cbind(msey, msey.unshrunk, msey.oracular)
colnames(msey2.table) = c("GWL", "GWL-U", "Oracle")
msey2bold = matrix(FALSE, nrow=dim(msey2.table)[1], ncol=dim(msey2.table)[2])
for (i in 1:(dim(msey2.table)[1])) {msey2bold[i,order(msey2.table[i,])[1]] = TRUE}
msey2ital = matrix(FALSE, nrow=dim(msey2.table)[1], ncol=dim(msey2.table)[2])
for (i in 1:(dim(msey2.table)[1])) {msey2ital[i,order(msey2.table[i,])[2]] = TRUE}
xtable.printbold(xtable(msey2.table, digits=3, caption="Mean squared error of estimates for the response variable $y$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEY}"), which.bold=msey2bold, which.ital=msey2ital, include.rownames=FALSE)

msey.table = cbind(msey, msey.unshrunk, msey.precon, msey.unshrunk.precon, msey.oracular)
colnames(msey.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
mseybold = matrix(FALSE, nrow=dim(msey.table)[1], ncol=dim(msey.table)[2])
for (i in 1:(dim(msey.table)[1])) {mseybold[i,order(msey.table[i,])[1]] = TRUE}
mseyital = matrix(FALSE, nrow=dim(msey.table)[1], ncol=dim(msey.table)[2])
for (i in 1:(dim(msey.table)[1])) {mseyital[i,order(msey.table[i,])[2]] = TRUE}
xtable.printbold(xtable(msey.table, digits=3), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE)

msey.loc.table = list()
for (l in 1:length(locs)) {
    msey.loc.table[[l]] = cbind(msey.loc[[l]][['GWL']], msey.loc[[l]][['unshrunk']], msey.loc[[l]][['precon']], msey.loc[[l]][['unshrunk-precon']], msey.loc[[l]][['oracular']])
    colnames(msey.loc.table[[l]]) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
    mseybold = matrix(FALSE, nrow=dim(msey.loc.table[[l]])[1], ncol=dim(msey.loc.table[[l]])[2])
    for (i in 1:(dim(msey.loc.table[[l]])[1])) {mseybold[i,order(msey.loc.table[[l]][i,])[1]] = TRUE}
    mseyital = matrix(FALSE, nrow=dim(msey.loc.table[[l]])[1], ncol=dim(msey.loc.table[[l]])[2])
    for (i in 1:(dim(msey.loc.table[[l]])[1])) {mseyital[i,order(msey.loc.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(msey.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c'), caption=paste("Mean squared error of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}", sep="")), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE, hline.after=c(0))
}




locs = c(30, 318, 450, 613, 871)
bx.loc = list()
for (l in 1:length(locs)) {    
    bx.loc[[l]] = list()
    bx.loc[[l]][['GWL']] = vector()
    bx.loc[[l]][['oracular']] = vector()
    bx.loc[[l]][['precon']] = vector()
    bx.loc[[l]][['unshrunk']] = vector()
    bx.loc[[l]][['unshrunk-precon']] = vector()
}

b2 = vector()
b2.oracular = vector()
b2.precon = vector()
b2.unshrunk = vector()
b2.unshrunk.precon = vector()

b2.nonzero = vector()
b2.oracular.nonzero = vector()
b2.precon.nonzero = vector()
b2.unshrunk.nonzero = vector()
b2.unshrunk.precon.nonzero = vector()
for (s in 1:18) {
	b2 = c(b2, mean(sapply(Map(function(z) {mean(sapply(X1.err[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.oracular = c(b2.oracular, mean(sapply(Map(function(z) {mean(sapply(X1.err.oracular[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.precon = c(b2.precon, mean(sapply(Map(function(z) {mean(sapply(X1.err.precon[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.unshrunk = c(b2.unshrunk, mean(sapply(Map(function(z) {mean(sapply(X1.err.unshrunk[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.unshrunk.precon = c(b2.unshrunk.precon, mean(sapply(Map(function(z) {mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[z]}))**2}, 1:900), identity)))

	b2.nonzero = c(b2.nonzero, mean(sapply(Map(function(z) {mean(sapply(X1.err[[s]], function(x) {x[z]}))**2}, which(B[['X1']]!=0)), identity)))
	b2.oracular.nonzero = c(b2.oracular.nonzero, mean(sapply(Map(function(z) {mean(sapply(X1.err.oracular[[s]], function(x) {x[z]}))**2}, which(B[['X1']]!=0)), identity)))
	b2.precon.nonzero = c(b2.precon.nonzero, mean(sapply(Map(function(z) {mean(sapply(X1.err.precon[[s]], function(x) {x[z]}))**2}, which(B[['X1']]!=0)), identity)))
	b2.unshrunk.nonzero = c(b2.unshrunk.nonzero, mean(sapply(Map(function(z) {mean(sapply(X1.err.unshrunk[[s]], function(x) {x[z]}))**2}, which(B[['X1']]!=0)), identity)))
	b2.unshrunk.precon.nonzero = c(b2.unshrunk.precon.nonzero, mean(sapply(Map(function(z) {mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[z]}))**2}, which(B[['X1']]!=0)), identity)))

    for (l in 1:length(locs)) {
        bx.loc[[l]][['GWL']] = c(bx.loc[[l]][['GWL']], mean(sapply(X1.err[[s]], function(x) {x[locs[l]]})))
        bx.loc[[l]][['oracular']] = c(bx.loc[[l]][['oracular']], mean(sapply(X1.err.oracular[[s]], function(x) {x[locs[l]]})))
        bx.loc[[l]][['precon']] = c(bx.loc[[l]][['precon']], mean(sapply(X1.err.precon[[s]], function(x) {x[locs[l]]})))
        bx.loc[[l]][['unshrunk']] = c(bx.loc[[l]][['unshrunk']], mean(sapply(X1.err.unshrunk[[s]], function(x) {x[locs[l]]})))
        bx.loc[[l]][['unshrunk-precon']] = c(bx.loc[[l]][['unshrunk-precon']], mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[locs[l]]})))
    }
}

b22.table = cbind(b2, b2.unshrunk, b2.oracular)
colnames(b22.table) = c("GWL", "GWL-U", "Oracle")
b22bold = matrix(FALSE, nrow=dim(b22.table)[1], ncol=dim(b22.table)[2])
for (i in 1:(dim(b22.table)[1])) {b22bold[i,order(b22.table[i,])[1]] = TRUE}
b22ital = matrix(FALSE, nrow=dim(b22.table)[1], ncol=dim(b22.table)[2])
for (i in 1:(dim(b22.table)[1])) {b22ital[i,order(b22.table[i,])[2]] = TRUE}
xtable.printbold(xtable(b22.table, digits=4, caption="Squared bias of estimates for $\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{bias}"), which.bold=b22bold, which.ital=b22ital, include.rownames=FALSE)

b22.nonzero.table = cbind(b2.nonzero, b2.unshrunk.nonzero, b2.oracular.nonzero)
colnames(b22.nonzero.table) = c("GWL", "GWL-U", "Oracle")
b22bold.nonzero = matrix(FALSE, nrow=dim(b22.nonzero.table)[1], ncol=dim(b22.nonzero.table)[2])
for (i in 1:(dim(b22.nonzero.table)[1])) {b22bold.nonzero[i,order(b22.nonzero.table[i,])[1]] = TRUE}
b22ital.nonzero = matrix(FALSE, nrow=dim(b22.nonzero.table)[1], ncol=dim(b22.nonzero.table)[2])
for (i in 1:(dim(b22.nonzero.table)[1])) {b22ital.nonzero[i,order(b22.nonzero.table[i,])[2]] = TRUE}
xtable.printbold(xtable(b22.nonzero.table, digits=4, caption="Squared bias of estimates for $\beta_1$ at locations where $\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{bias-nonzero}"), which.bold=b22bold.nonzero, which.ital=b22ital.nonzero, include.rownames=FALSE)


b2.table = cbind(b2, b2.unshrunk, b2.precon, b2.unshrunk.precon, b2.oracular)
colnames(b2.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
b2bold = matrix(FALSE, nrow=dim(b2.table)[1], ncol=dim(b2.table)[2])
for (i in 1:(dim(b2.table)[1])) {b2bold[i,order(b2.table[i,])[1]] = TRUE}
b2ital = matrix(FALSE, nrow=dim(b2.table)[1], ncol=dim(b2.table)[2])
for (i in 1:(dim(b2.table)[1])) {b2ital[i,order(b2.table[i,])[2]] = TRUE}
xtable.printbold(xtable(b2.table, digits=4), which.bold=b2bold, which.ital=b2ital, include.rownames=FALSE)

b2.nonzero.table = cbind(b2.nonzero, b2.unshrunk.nonzero, b2.precon.nonzero, b2.unshrunk.precon.nonzero, b2.oracular.nonzero)
colnames(b2.nonzero.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
b2bold.nonzero = matrix(FALSE, nrow=dim(b2.nonzero.table)[1], ncol=dim(b2.nonzero.table)[2])
for (i in 1:(dim(b2.nonzero.table)[1])) {b2bold.nonzero[i,order(b2.nonzero.table[i,])[1]] = TRUE}
b2ital.nonzero = matrix(FALSE, nrow=dim(b2.nonzero.table)[1], ncol=dim(b2.nonzero.table)[2])
for (i in 1:(dim(b2.nonzero.table)[1])) {b2ital.nonzero[i,order(b2.nonzero.table[i,])[2]] = TRUE}
xtable.printbold(xtable(b2.nonzero.table, digits=4), which.bold=b2bold.nonzero, which.ital=b2ital.nonzero, include.rownames=FALSE, hline.after=c(0))

bx.loc.table = list()
for (l in 1:length(locs)) {
    bx.loc.table[[l]] = cbind(bx.loc[[l]][['GWL']], bx.loc[[l]][['unshrunk']], bx.loc[[l]][['precon']], bx.loc[[l]][['unshrunk-precon']], bx.loc[[l]][['oracular']])
    colnames(bx.loc.table[[l]]) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
    bxbold = matrix(FALSE, nrow=dim(bx.loc.table[[l]])[1], ncol=dim(bx.loc.table[[l]])[2])
    for (i in 1:(dim(bx.loc.table[[l]])[1])) {bxbold[i,order(abs(bx.loc.table[[l]][i,]))[1]] = TRUE}
    bxital = matrix(FALSE, nrow=dim(bx.loc.table[[l]])[1], ncol=dim(bx.loc.table[[l]])[2])
    for (i in 1:(dim(bx.loc.table[[l]])[1])) {bxital[i,order(abs(bx.loc.table[[l]][i,]))[2]] = TRUE}
    xtable.printbold(xtable(bx.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c'), caption=paste("Bias of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}", sep="")), which.bold=bxbold, which.ital=bxital, include.rownames=FALSE, hline.after=c(0))
}






locs = c(30, 318, 450, 613, 871)
by.loc = list()
for (l in 1:length(locs)) {    
    by.loc[[l]] = list()
    by.loc[[l]][['GWL']] = vector()
    by.loc[[l]][['oracular']] = vector()
    by.loc[[l]][['precon']] = vector()
    by.loc[[l]][['unshrunk']] = vector()
    by.loc[[l]][['unshrunk-precon']] = vector()
}

by2 = vector()
by2.oracular = vector()
by2.precon = vector()
by2.unshrunk = vector()
by2.unshrunk.precon = vector()
for (s in 1:18) {
	by2 = c(by2, mean(sapply(Map(function(z) {mean(sapply(Y.err[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	by2.oracular = c(by2.oracular, mean(sapply(Map(function(z) {mean(sapply(Y.err.oracular[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	by2.precon = c(by2.precon, mean(sapply(Map(function(z) {mean(sapply(Y.err.precon[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	by2.unshrunk = c(by2.unshrunk, mean(sapply(Map(function(z) {mean(sapply(Y.err.unshrunk[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	by2.unshrunk.precon = c(by2.unshrunk.precon, mean(sapply(Map(function(z) {mean(sapply(Y.err.unshrunk.precon[[s]], function(x) {x[z]}))**2}, 1:900), identity)))

    for (l in 1:length(locs)) {
        by.loc[[l]][['GWL']] = c(by.loc[[l]][['GWL']], mean(sapply(Y.err[[s]], function(x) {x[locs[l]]})))
        by.loc[[l]][['oracular']] = c(by.loc[[l]][['oracular']], mean(sapply(Y.err.oracular[[s]], function(x) {x[locs[l]]})))
        by.loc[[l]][['precon']] = c(by.loc[[l]][['precon']], mean(sapply(Y.err.precon[[s]], function(x) {x[locs[l]]})))
        by.loc[[l]][['unshrunk']] = c(by.loc[[l]][['unshrunk']], mean(sapply(Y.err.unshrunk[[s]], function(x) {x[locs[l]]})))
        by.loc[[l]][['unshrunk-precon']] = c(by.loc[[l]][['unshrunk-precon']], mean(sapply(Y.err.unshrunk.precon[[s]], function(x) {x[locs[l]]})))
    }
}
by22.table = cbind(by2, by2.unshrunk, by2.oracular)
colnames(by22.table) = c("GWL", "GWL-U", "Oracle")
by22bold = matrix(FALSE, nrow=dim(by22.table)[1], ncol=dim(by22.table)[2])
for (i in 1:(dim(by22.table)[1])) {by22bold[i,order(by22.table[i,])[1]] = TRUE}
by22ital = matrix(FALSE, nrow=dim(by22.table)[1], ncol=dim(by22.table)[2])
for (i in 1:(dim(by22.table)[1])) {by22ital[i,order(by22.table[i,])[2]] = TRUE}
xtable.printbold(xtable(by22.table, digits=4, caption="Squared bias of estimates for the response variable $y$ (\\textbf{minimum}, \\emph{next best}).\\label{biasy}"), which.bold=by22bold, which.ital=by22ital, include.rownames=FALSE)


by2.table = cbind(by2, by2.unshrunk, by2.precon, by2.unshrunk.precon, by2.oracular)
colnames(by2.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
by2bold = matrix(FALSE, nrow=dim(by2.table)[1], ncol=dim(by2.table)[2])
for (i in 1:(dim(by2.table)[1])) {by2bold[i,order(by2.table[i,])[1]] = TRUE}
by2ital = matrix(FALSE, nrow=dim(by2.table)[1], ncol=dim(by2.table)[2])
for (i in 1:(dim(by2.table)[1])) {by2ital[i,order(by2.table[i,])[2]] = TRUE}
xtable.printbold(xtable(by2.table, digits=4), which.bold=by2bold, which.ital=by2ital, include.rownames=FALSE)

by.loc.table = list()
for (l in 1:length(locs)) {
    by.loc.table[[l]] = cbind(by.loc[[l]][['GWL']], by.loc[[l]][['unshrunk']], by.loc[[l]][['precon']], by.loc[[l]][['unshrunk-precon']], by.loc[[l]][['oracular']])
    colnames(by.loc.table[[l]]) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
    bybold = matrix(FALSE, nrow=dim(by.loc.table[[l]])[1], ncol=dim(by.loc.table[[l]])[2])
    for (i in 1:(dim(by.loc.table[[l]])[1])) {bybold[i,order(abs(by.loc.table[[l]][i,]))[1]] = TRUE}
    byital = matrix(FALSE, nrow=dim(by.loc.table[[l]])[1], ncol=dim(by.loc.table[[l]])[2])
    for (i in 1:(dim(by.loc.table[[l]])[1])) {byital[i,order(abs(by.loc.table[[l]][i,]))[2]] = TRUE}
    xtable.printbold(xtable(by.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c'), caption=paste("Bias of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).", sep="")), which.bold=bybold, which.ital=byital, include.rownames=FALSE, hline.after=c(0))
}




locs = c(30, 318, 450, 613, 871)
varx.loc = list()
for (l in 1:length(locs)) {    
    varx.loc[[l]] = list()
    varx.loc[[l]][['GWL']] = vector()
    varx.loc[[l]][['oracular']] = vector()
    varx.loc[[l]][['precon']] = vector()
    varx.loc[[l]][['unshrunk']] = vector()
    varx.loc[[l]][['unshrunk-precon']] = vector()
}

varx = vector()
varx.oracular = vector()
varx.precon = vector()
varx.unshrunk = vector()
varx.unshrunk.precon = vector()

varx.nonzero = vector()
varx.oracular.nonzero = vector()
varx.precon.nonzero = vector()
varx.unshrunk.nonzero = vector()
varx.unshrunk.precon.nonzero = vector()
for (s in 1:18) {
	varx = c(varx, mean(sapply(Map(function(z) {var(sapply(X1.err[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.oracular = c(varx.oracular, mean(sapply(Map(function(z) {var(sapply(X1.err.oracular[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.precon = c(varx.precon, mean(sapply(Map(function(z) {var(sapply(X1.err.precon[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.unshrunk = c(varx.unshrunk, mean(sapply(Map(function(z) {var(sapply(X1.err.unshrunk[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.unshrunk.precon = c(varx.unshrunk.precon, mean(sapply(Map(function(z) {var(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[z]}))}, 1:900), identity)))

	varx.nonzero = c(varx.nonzero, mean(sapply(Map(function(z) {var(sapply(X1.err[[s]], function(x) {x[z]}))}, which(B[['X1']]!=0)), identity)))
	varx.oracular.nonzero = c(varx.oracular.nonzero, mean(sapply(Map(function(z) {var(sapply(X1.err.oracular[[s]], function(x) {x[z]}))}, which(B[['X1']]!=0)), identity)))
	varx.precon.nonzero = c(varx.precon.nonzero, mean(sapply(Map(function(z) {var(sapply(X1.err.precon[[s]], function(x) {x[z]}))}, which(B[['X1']]!=0)), identity)))
	varx.unshrunk.nonzero = c(varx.unshrunk.nonzero, mean(sapply(Map(function(z) {var(sapply(X1.err.unshrunk[[s]], function(x) {x[z]}))}, which(B[['X1']]!=0)), identity)))
	varx.unshrunk.precon.nonzero = c(varx.unshrunk.precon.nonzero, mean(sapply(Map(function(z) {var(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[z]}))}, which(B[['X1']]!=0)), identity)))

    for (l in 1:length(locs)) {
        varx.loc[[l]][['GWL']] = c(varx.loc[[l]][['GWL']], var(sapply(X1.err[[s]], function(x) {x[locs[l]]})))
        varx.loc[[l]][['oracular']] = c(varx.loc[[l]][['oracular']], var(sapply(X1.err.oracular[[s]], function(x) {x[locs[l]]})))
        varx.loc[[l]][['precon']] = c(varx.loc[[l]][['precon']], var(sapply(X1.err.precon[[s]], function(x) {x[locs[l]]})))
        varx.loc[[l]][['unshrunk']] = c(varx.loc[[l]][['unshrunk']], var(sapply(X1.err.unshrunk[[s]], function(x) {x[locs[l]]})))
        varx.loc[[l]][['unshrunk-precon']] = c(varx.loc[[l]][['unshrunk-precon']], var(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[locs[l]]})))
    }
}
varx2.table = cbind(varx, varx.unshrunk, varx.oracular)
colnames(varx2.table) = c("GWL", "GWL-U", "Oracle")
varx2bold = matrix(FALSE, nrow=dim(varx2.table)[1], ncol=dim(varx2.table)[2])
for (i in 1:(dim(varx2.table)[1])) {varx2bold[i,order(varx2.table[i,])[1]] = TRUE}
varx2ital = matrix(FALSE, nrow=dim(varx2.table)[1], ncol=dim(varx2.table)[2])
for (i in 1:(dim(varx2.table)[1])) {varx2ital[i,order(varx2.table[i,])[2]] = TRUE}
xtable.printbold(xtable(varx2.table, digits=4, caption="Variance of estimates for $\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{varx}"), which.bold=varx2bold, which.ital=varx2ital, include.rownames=FALSE)

varx2.nonzero.table = cbind(varx.nonzero, varx.unshrunk.nonzero, varx.oracular.nonzero)
colnames(varx2.nonzero.table) = c("GWL", "GWL-U", "Oracle")
varx2bold.nonzero = matrix(FALSE, nrow=dim(varx2.nonzero.table)[1], ncol=dim(varx2.nonzero.table)[2])
for (i in 1:(dim(varx2.nonzero.table)[1])) {varx2bold.nonzero[i,order(varx2.nonzero.table[i,])[1]] = TRUE}
varx2ital.nonzero = matrix(FALSE, nrow=dim(varx2.nonzero.table)[1], ncol=dim(varx2.nonzero.table)[2])
for (i in 1:(dim(varx2.nonzero.table)[1])) {varx2ital.nonzero[i,order(varx2.nonzero.table[i,])[2]] = TRUE}
xtable.printbold(xtable(varx2.nonzero.table, digits=4, caption="Variance of estimates for $\beta_1$ at locations where $\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{varx-nonzero}"), which.bold=varx2bold.nonzero, which.ital=varx2ital.nonzero, include.rownames=FALSE)


varx.table = cbind(varx, varx.unshrunk, varx.precon, varx.unshrunk.precon, varx.oracular)
colnames(varx.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
varxbold = matrix(FALSE, nrow=dim(varx.table)[1], ncol=dim(varx.table)[2])
for (i in 1:(dim(varx.table)[1])) {varxbold[i,order(varx.table[i,])[1]] = TRUE}
varxital = matrix(FALSE, nrow=dim(varx.table)[1], ncol=dim(varx.table)[2])
for (i in 1:(dim(varx.table)[1])) {varxital[i,order(varx.table[i,])[2]] = TRUE}
xtable.printbold(xtable(varx.table, digits=4), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE)

varx.nonzero.table = cbind(varx.nonzero, varx.unshrunk.nonzero, varx.precon.nonzero, varx.unshrunk.precon.nonzero, varx.oracular.nonzero)
colnames(varx.nonzero.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
varxbold.nonzero = matrix(FALSE, nrow=dim(varx.nonzero.table)[1], ncol=dim(varx.nonzero.table)[2])
for (i in 1:(dim(varx.nonzero.table)[1])) {varxbold.nonzero[i,order(varx.nonzero.table[i,])[1]] = TRUE}
varxital.nonzero = matrix(FALSE, nrow=dim(varx.nonzero.table)[1], ncol=dim(varx.nonzero.table)[2])
for (i in 1:(dim(varx.nonzero.table)[1])) {varxital.nonzero[i,order(varx.nonzero.table[i,])[2]] = TRUE}
xtable.printbold(xtable(varx.nonzero.table, digits=4), which.bold=varxbold.nonzero, which.ital=varxital.nonzero, include.rownames=FALSE, hline.after=c(0))

varx.loc.table = list()
for (l in 1:length(locs)) {
    varx.loc.table[[l]] = cbind(varx.loc[[l]][['GWL']], varx.loc[[l]][['unshrunk']], varx.loc[[l]][['precon']], varx.loc[[l]][['unshrunk-precon']], varx.loc[[l]][['oracular']])
    colnames(varx.loc.table[[l]]) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
    varxbold = matrix(FALSE, nrow=dim(varx.loc.table[[l]])[1], ncol=dim(varx.loc.table[[l]])[2])
    for (i in 1:(dim(varx.loc.table[[l]])[1])) {varxbold[i,order(varx.loc.table[[l]][i,])[1]] = TRUE}
    varxital = matrix(FALSE, nrow=dim(varx.loc.table[[l]])[1], ncol=dim(varx.loc.table[[l]])[2])
    for (i in 1:(dim(varx.loc.table[[l]])[1])) {varxital[i,order(varx.loc.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(varx.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c'), caption=paste("Variance of estimates for $\\beta_1$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).", sep="")), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE, hline.after=c(0))
}




locs = c(30, 318, 450, 613, 871)
vary.loc = list()
for (l in 1:length(locs)) {    
    vary.loc[[l]] = list()
    vary.loc[[l]][['GWL']] = vector()
    vary.loc[[l]][['oracular']] = vector()
    vary.loc[[l]][['precon']] = vector()
    vary.loc[[l]][['unshrunk']] = vector()
    vary.loc[[l]][['unshrunk-precon']] = vector()
}

vary = vector()
vary.oracular = vector()
vary.precon = vector()
vary.unshrunk = vector()
vary.unshrunk.precon = vector()
for (s in 1:18) {
	vary = c(vary, mean(sapply(Map(function(z) {var(sapply(Y.err[[s]], function(x) {x[z]}))}, 1:900), identity)))
	vary.oracular = c(vary.oracular, mean(sapply(Map(function(z) {var(sapply(Y.err.oracular[[s]], function(x) {x[z]}))}, 1:900), identity)))
	vary.precon = c(vary.precon, mean(sapply(Map(function(z) {var(sapply(Y.err.precon[[s]], function(x) {x[z]}))}, 1:900), identity)))
	vary.unshrunk = c(vary.unshrunk, mean(sapply(Map(function(z) {var(sapply(Y.err.unshrunk[[s]], function(x) {x[z]}))}, 1:900), identity)))
	vary.unshrunk.precon = c(vary.unshrunk.precon, mean(sapply(Map(function(z) {var(sapply(Y.err.unshrunk.precon[[s]], function(x) {x[z]}))}, 1:900), identity)))

    for (l in 1:length(locs)) {
        vary.loc[[l]][['GWL']] = c(vary.loc[[l]][['GWL']], var(sapply(Y.err[[s]], function(x) {x[locs[l]]})))
        vary.loc[[l]][['oracular']] = c(vary.loc[[l]][['oracular']], var(sapply(Y.err.oracular[[s]], function(x) {x[locs[l]]})))
        vary.loc[[l]][['precon']] = c(vary.loc[[l]][['precon']], var(sapply(Y.err.precon[[s]], function(x) {x[locs[l]]})))
        vary.loc[[l]][['unshrunk']] = c(vary.loc[[l]][['unshrunk']], var(sapply(Y.err.unshrunk[[s]], function(x) {x[locs[l]]})))
        vary.loc[[l]][['unshrunk-precon']] = c(vary.loc[[l]][['unshrunk-precon']], var(sapply(Y.err.unshrunk.precon[[s]], function(x) {x[locs[l]]})))
    }
}
vary2.table = cbind(vary, vary.unshrunk, vary.oracular)
colnames(vary2.table) = c("GWL", "GWL-U", "Oracle")
vary2bold = matrix(FALSE, nrow=dim(vary2.table)[1], ncol=dim(vary2.table)[2])
for (i in 1:(dim(vary2.table)[1])) {vary2bold[i,order(vary2.table[i,])[1]] = TRUE}
vary2ital = matrix(FALSE, nrow=dim(vary2.table)[1], ncol=dim(vary2.table)[2])
for (i in 1:(dim(vary2.table)[1])) {vary2ital[i,order(vary2.table[i,])[2]] = TRUE}
xtable.printbold(xtable(vary2.table, digits=4, caption="Variance of estimates for the response variable $y$ (\\textbf{minimum}, \\emph{next best}).\\label{vary}"), which.bold=vary2bold, which.ital=vary2ital, include.rownames=FALSE)

vary.table = cbind(vary, vary.unshrunk, vary.precon, vary.unshrunk.precon, vary.oracular)
colnames(vary.table) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
varybold = matrix(FALSE, nrow=dim(vary.table)[1], ncol=dim(vary.table)[2])
for (i in 1:(dim(vary.table)[1])) {varybold[i,order(vary.table[i,])[1]] = TRUE}
varyital = matrix(FALSE, nrow=dim(vary.table)[1], ncol=dim(vary.table)[2])
for (i in 1:(dim(vary.table)[1])) {varyital[i,order(vary.table[i,])[2]] = TRUE}
xtable.printbold(xtable(vary.table, digits=4), which.bold=varybold, which.ital=varyital, include.rownames=FALSE)

vary.loc.table = list()
for (l in 1:length(locs)) {
    vary.loc.table[[l]] = cbind(vary.loc[[l]][['GWL']], vary.loc[[l]][['unshrunk']], vary.loc[[l]][['precon']], vary.loc[[l]][['unshrunk-precon']], vary.loc[[l]][['oracular']])
    colnames(vary.loc.table[[l]]) = c("GWL", "GWL-U", "GWL-P", "GWL-P-U", "Oracle")
    varybold = matrix(FALSE, nrow=dim(vary.loc.table[[l]])[1], ncol=dim(vary.loc.table[[l]])[2])
    for (i in 1:(dim(vary.loc.table[[l]])[1])) {varybold[i,order(vary.loc.table[[l]][i,])[1]] = TRUE}
    varyital = matrix(FALSE, nrow=dim(vary.loc.table[[l]])[1], ncol=dim(vary.loc.table[[l]])[2])
    for (i in 1:(dim(vary.loc.table[[l]])[1])) {varyital[i,order(vary.loc.table[[l]][i,])[2]] = TRUE}
    xtable.printbold(xtable(vary.loc.table[[l]], digits=3, align=c('c','c','c','c','c','c'), caption=paste("Variance of estimates for $Y$ at location ", l, " (\\textbf{minimum}, \\emph{next best}).", sep="")), which.bold=varybold, which.ital=varyital, include.rownames=FALSE, hline.after=c(0))
}




vv = c('X1', 'X2', 'X3', 'X4', 'X5')
errs.zeros.denom = rep(0, N**2)
errs.nonzeros.denom = rep(0, N**2)
for (v in vv) {
    errs.zeros.denom = errs.zeros.denom + ifelse(B[[v]]==0, 1, 0)
    errs.nonzeros.denom = errs.nonzeros.denom + ifelse(B[[v]]!=0, 1, 0)
}

errs.zeros = list()
errs.nonzeros = list()
for (s in 1:18) {
    err.z = rep(0, N**2)
    err.nz = rep(0, N**2)

    for (v in vv) {
	    err.z = err.z + apply( selection[[s]][[v]], 2, function(x) {ifelse(B[[v]]==0, x, 0)})
        err.nz = err.nz + apply( selection[[s]][[v]], 2, function(x) {ifelse(B[[v]]!=0, 1-x, 0)})
    }

	errs.zeros[[s]] = rowSums(err.z)
    errs.nonzeros[[s]] = rowSums(err.nz)
}

errs.zeros.precon = list()
errs.nonzeros.precon = list()
for (s in 1:18) {
    err.z = rep(0, N**2)
    err.nz = rep(0, N**2)

    for (v in vv) {
	    err.z = err.z + apply( selection.precon[[s]][[v]], 2, function(x) {ifelse(B[[v]]==0, x, 0)})
        err.nz = err.nz + apply( selection.precon[[s]][[v]], 2, function(x) {ifelse(B[[v]]!=0, 1-x, 0)})
    }

	errs.zeros.precon[[s]] = rowSums(err.z)
    errs.nonzeros.precon[[s]] = rowSums(err.nz)
}


perfect.selection = data.frame()
for (s in 1:18) {
	perfect.selection = rbind(perfect.selection, c(sum(errs.zeros[[s]]/100)/sum(errs.zeros.denom), sum(errs.nonzeros[[s]]/100) / sum(errs.nonzeros.denom), sum(errs.zeros.precon[[s]]/100) / sum(errs.zeros.denom), sum(errs.nonzeros.precon[[s]]/100) / sum(errs.nonzeros.denom)))
}
colnames(perfect.selection) = c("original", "preconditioned")

selectbold = matrix(FALSE, nrow=dim(perfect.selection)[1], ncol=dim(perfect.selection)[2])
for (i in 1:(dim(perfect.selection)[1])) {selectbold[i,order(perfect.selection[i,], decreasing=TRUE)[1]] = TRUE}

xtable.printbold(xtable(perfect.selection, digits=3), which.bold=selectbold, include.rownames=FALSE)





tau = rep(c(0.03, 0.1), each=9)
rho = rep(rep(c(0, 0.5, 0.8), each=3), times=2)
sigma.tau = rep(c(0, 0.03, 0.1), times=6)

t.x = list()
t.x[[1]] = 1:9
t.x[[2]] = 10:18

rho.2 = list()
rho.2[[1]] = c(1:3, 10:12)
rho.2[[2]] = c(4:6, 13:15)
rho.2[[3]] = c(7:9, 16:18)

t.e = list()
t.e[[1]] = c(1,4,7,10,13,16)
t.e[[2]] = c(2,5,8,11,14,17)
t.e[[3]] = c(3,6,9,12,15,18)


selection.table = list()
selection.table[['t.x']] = c(mean(perfect.selection[t.x[[1]],1]), mean(perfect.selection[t.x[[2]],1]))
selection.table[['t.e']] = c(mean(perfect.selection[t.e[[1]],1]), mean(perfect.selection[t.e[[2]],1]), mean(perfect.selection[t.e[[3]],1]))
selection.table[['rho']] = c(mean(perfect.selection[rho.2[[1]],1]), mean(perfect.selection[rho.2[[2]],1]), mean(perfect.selection[rho.2[[3]],1]))

tau.table = data.frame(unique(tau), selection.table[['t.x']])
tau.sigma.table = data.frame(unique(sigma.tau), selection.table[['t.e']])
rho.table = data.frame(unique(rho), selection.table[['rho']])

colnames(tau.table) = c("$\\tau_x$", "Frequency")
colnames(tau.sigma.table) = c("$\\tau_\\sigma$", "Frequency")
colnames(rho.table) = c("$\\rho$", "Frequency")

print.xtable(xtable(tau.table, digits=c(0,2,2), caption="Frequency of perfect selection for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), sanitize.colnames.function=function(x){x})
print.xtable(xtable(tau.sigma.table, digits=c(0,2,2), caption="Frequency of perfect selection for different settings of $\\tau_{\\sigma}$."), include.rownames=FALSE, hline.after=c(0), sanitize.colnames.function=function(x){x})
print.xtable(xtable(rho.table, digits=c(0,2,2), caption="Frequency of perfect selection for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), sanitize.colnames.function=function(x){x})


msex.table = list()
msex.table[['t.x']] = cbind(unique(tau), rbind(apply(mse2[t.x[[1]],],2,mean), apply(mse2[t.x[[2]],],2,mean)))
msex.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(mse2[t.e[[1]],],2,mean), apply(mse2[t.e[[2]],],2,mean), apply(mse2[t.e[[3]],],2,mean)))
msex.table[['rho']] = cbind(unique(rho), rbind(apply(mse2[rho.2[[1]],],2,mean), apply(mse2[rho.2[[2]],],2,mean), apply(mse2[rho.2[[3]],],2,mean)))
colnames(msex.table[['t.x']]) = c("$\\tau_x$", "GWL", "GWL-U", "Oracle")
colnames(msex.table[['t.e']]) = c("$\\tau_\\sigma$", "GWL", "GWL-U", "Oracle")
colnames(msex.table[['rho']]) = c("$\\rho$", "GWL", "GWL-U", "Oracle")

bold.msex = list()
bold.msex[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.msex[['t.x']][i,1+which.min(msex.table[['t.x']][i,2:4])] = TRUE}
bold.msex[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['t.e']][i,1+which.min(msex.table[['t.e']][i,2:4])] = TRUE}
bold.msex[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['rho']][i,1+which.min(msex.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(msex.table[['t.x']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.x']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(msex.table[['t.e']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.e']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(msex.table[['rho']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['rho']], sanitize.colnames.function=function(x){x})

msex.table.nz = list()
msex.table.nz[['t.x']] = cbind(unique(tau), rbind(apply(mse2.nonzero[t.x[[1]],],2,mean), apply(mse2.nonzero[t.x[[2]],],2,mean)))
msex.table.nz[['t.e']] = cbind(unique(sigma.tau), rbind(apply(mse2.nonzero[t.e[[1]],],2,mean), apply(mse2.nonzero[t.e[[2]],],2,mean), apply(mse2.nonzero[t.e[[3]],],2,mean)))
msex.table.nz[['rho']] = cbind(unique(rho), rbind(apply(mse2.nonzero[rho.2[[1]],],2,mean), apply(mse2.nonzero[rho.2[[2]],],2,mean), apply(mse2.nonzero[rho.2[[3]],],2,mean)))
colnames(msex.table.nz[['t.x']]) = c("$\\tau_x$", "GWL", "GWL-U", "Oracle")
colnames(msex.table.nz[['t.e']]) = c("$\\tau_\\sigma$", "GWL", "GWL-U", "Oracle")
colnames(msex.table.nz[['rho']]) = c("$\\rho$", "GWL", "GWL-U", "Oracle")

bold.msex.nz = list()
bold.msex.nz[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.msex.nz[['t.x']][i,1+which.min(msex.table[['t.x']][i,2:4])] = TRUE}
bold.msex[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['t.e']][i,1+which.min(msex.table[['t.e']][i,2:4])] = TRUE}
bold.msex[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['rho']][i,1+which.min(msex.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(msex.table[['t.x']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.x']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(msex.table[['t.e']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.e']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(msex.table[['rho']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['rho']], sanitize.colnames.function=function(x){x})



varx.table = list()
varx.table[['t.x']] = cbind(unique(tau), rbind(apply(varx2[t.x[[1]],],2,mean), apply(varx2[t.x[[2]],],2,mean)))
varx.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(varx2[t.e[[1]],],2,mean), apply(varx2[t.e[[2]],],2,mean), apply(varx2[t.e[[3]],],2,mean)))
varx.table[['rho']] = cbind(unique(rho), rbind(apply(varx2[rho.2[[1]],],2,mean), apply(varx2[rho.2[[2]],],2,mean), apply(varx2[rho.2[[3]],],2,mean)))
colnames(varx.table[['t.x']]) = c("$\\tau_x$", "GWL", "GWL-U", "Oracle")
colnames(varx.table[['t.e']]) = c("$\\tau_\\sigma$", "GWL", "GWL-U", "Oracle")
colnames(varx.table[['rho']]) = c("$\\rho$", "GWL", "GWL-U", "Oracle")

bold.varx = list()
bold.varx[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.varx[['t.x']][i,1+which.min(varx.table[['t.x']][i,2:4])] = TRUE}
bold.varx[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.varx[['t.e']][i,1+which.min(varx.table[['t.e']][i,2:4])] = TRUE}
bold.varx[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.varx[['rho']][i,1+which.min(varx.table[['rho']][i,2:4])] = TRUE}
xtable.printbold(xtable(varx.table[['t.x']], digits=c(0,2,3,3,3), caption="Variance of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.varx[['t.x']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(varx.table[['t.e']], digits=c(0,2,3,3,3), caption="Variance of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.varx[['t.e']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(varx.table[['rho']], digits=c(0,2,3,3,3), caption="Variance of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.varx[['rho']], sanitize.colnames.function=function(x){x})


bx.table = list()
bx.table[['t.x']] = cbind(unique(tau), rbind(apply(b22[t.x[[1]],],2,mean), apply(b22[t.x[[2]],],2,mean)))
bx.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(b22[t.e[[1]],],2,mean), apply(b22[t.e[[2]],],2,mean), apply(b22[t.e[[3]],],2,mean)))
bx.table[['rho']] = cbind(unique(rho), rbind(apply(b22[rho.2[[1]],],2,mean), apply(b22[rho.2[[2]],],2,mean), apply(b22[rho.2[[3]],],2,mean)))
colnames(bx.table[['t.x']]) = c("$\\tau_x$", "GWL", "GWL-U", "Oracle")
colnames(bx.table[['t.e']]) = c("$\\tau_\\sigma$", "GWL", "GWL-U", "Oracle")
colnames(bx.table[['rho']]) = c("$\\rho$", "GWL", "GWL-U", "Oracle")

bold.bx = list()
bold.bx[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.bx[['t.x']][i,1+which.min(bx.table[['t.x']][i,2:4])] = TRUE}
bold.bx[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.bx[['t.e']][i,1+which.min(bx.table[['t.e']][i,2:4])] = TRUE}
bold.bx[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.bx[['rho']][i,1+which.min(bx.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(bx.table[['t.x']], digits=c(0,2,3,3,3), caption="Squared bias of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.bx[['t.x']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(bx.table[['t.e']], digits=c(0,2,3,3,3), caption="Squared bias of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.bx[['t.e']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(bx.table[['rho']], digits=c(0,2,3,3,3), caption="Squared bias of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.bx[['rho']], sanitize.colnames.function=function(x){x})




by.table = list()
by.table[['t.x']] = cbind(unique(tau), rbind(apply(by22[t.x[[1]],],2,mean), apply(by22[t.x[[2]],],2,mean)))
by.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(by22[t.e[[1]],],2,mean), apply(by22[t.e[[2]],],2,mean), apply(by22[t.e[[3]],],2,mean)))
by.table[['rho']] = cbind(unique(rho), rbind(apply(by22[rho.2[[1]],],2,mean), apply(by22[rho.2[[2]],],2,mean), apply(by22[rho.2[[3]],],2,mean)))
colnames(by.table[['t.x']]) = c("$\\tau_x$", "GWL", "GWL-U", "Oracle")
colnames(by.table[['t.e']]) = c("$\\tau_\\sigma$", "GWL", "GWL-U", "Oracle")
colnames(by.table[['rho']]) = c("$\\rho$", "GWL", "GWL-U", "Oracle")

bold.by = list()
bold.by[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.by[['t.x']][i,1+which.min(by.table[['t.x']][i,2:4])] = TRUE}
bold.by[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.by[['t.e']][i,1+which.min(by.table[['t.e']][i,2:4])] = TRUE}
bold.by[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.by[['rho']][i,1+which.min(by.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(by.table[['t.x']], digits=c(0,2,3,3,3), caption="Squared bias of the estimated response variable $y$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.by[['t.x']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(by.table[['t.e']], digits=c(0,2,3,3,3), caption="Squared bias of the estimated response variable $y$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.by[['t.e']], sanitize.colnames.function=function(x){x})
xtable.printbold(xtable(by.table[['rho']], digits=c(0,2,3,3,3), caption="Squared bias of the estimated response variable $y$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.by[['rho']], sanitize.colnames.function=function(x){x})



