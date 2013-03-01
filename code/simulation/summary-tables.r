


mse = vector()
mse.oracular = vector()
mse.precon = vector()
mse.unshrunk = vector()
mse.unshrunk.precon = vector()

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

	mse.nonzero = c(mse.nonzero, mean(sapply(X1.err[[s]], function(x) {mean(x[B[['X1']]!=0]**2)}), na.rm=TRUE))
    mse.oracular.nonzero = c(mse.oracular.nonzero, mean(sapply(X1.err.oracular[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
	mse.precon.nonzero = c(mse.precon.nonzero, mean(sapply(X1.err.precon[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
	mse.unshrunk.nonzero = c(mse.unshrunk.nonzero, mean(sapply(X1.err.unshrunk[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
	mse.unshrunk.precon.nonzero = c(mse.unshrunk.precon.nonzero, mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {mean(x[B[['X1']]!=0]**2)})))
}
mse2 = cbind(mse, mse.unshrunk, mse.oracular)
colnames(mse2) = c("GWL", "GWL-U", "Oracle")
mse2bold = matrix(FALSE, nrow=dim(mse2)[1], ncol=dim(mse2)[2])
for (i in 1:(dim(mse2)[1])) {mse2bold[i,order(mse2[i,])[1]] = TRUE}
mse2ital = matrix(FALSE, nrow=dim(mse2)[1], ncol=dim(mse2)[2])
for (i in 1:(dim(mse2)[1])) {mse2ital[i,order(mse2[i,])[2]] = TRUE}
xtable.printbold(xtable(mse2, digits=3, align=c('c','c','c','c'), caption="Mean squared error of estimates for $\\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX}"), which.bold=mse2bold, which.ital=mse2ital, include.rownames=FALSE, hline.after=c(0))

mse2.nonzero = cbind(mse.nonzero, mse.unshrunk.nonzero, mse.oracular.nonzero)
colnames(mse2.nonzero) = c("GWL", "GWL-U", "Oracle")
mse2bold.nonzero = matrix(FALSE, nrow=dim(mse2.nonzero)[1], ncol=dim(mse2.nonzero)[2])
for (i in 1:(dim(mse2.nonzero)[1])) {mse2bold.nonzero[i,order(mse2.nonzero[i,])[1]] = TRUE}
mse2ital.nonzero = matrix(FALSE, nrow=dim(mse2.nonzero)[1], ncol=dim(mse2.nonzero)[2])
for (i in 1:(dim(mse2.nonzero)[1])) {mse2ital.nonzero[i,order(mse2.nonzero[i,])[2]] = TRUE}
xtable.printbold(xtable(mse2.nonzero, digits=3, align=c('c','c','c','c'), caption="Mean squared error of estimates for $\\beta_1$ at locations where $\\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEX-nonzero}"), which.bold=mse2bold.nonzero, which.ital=mse2ital.nonzero, include.rownames=FALSE, hline.after=c(0))


mse = cbind(mse, mse.unshrunk, mse.precon, mse.unshrunk.precon, mse.oracular)
colnames(mse) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
msebold = matrix(FALSE, nrow=dim(mse)[1], ncol=dim(mse)[2])
for (i in 1:(dim(mse)[1])) {msebold[i,order(mse[i,])[1]] = TRUE}
mseital = matrix(FALSE, nrow=dim(mse)[1], ncol=dim(mse)[2])
for (i in 1:(dim(mse)[1])) {mseital[i,order(mse[i,])[2]] = TRUE}
xtable.printbold(xtable(mse, digits=3, align=c('c','c','c','c','c','c')), which.bold=msebold, which.ital=mseital, include.rownames=FALSE, hline.after=c(0))

mse.nonzero = cbind(mse.nonzero, mse.unshrunk.nonzero, mse.precon.nonzero, mse.unshrunk.precon.nonzero, mse.oracular.nonzero)
colnames(mse.nonzero) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
msebold.nonzero = matrix(FALSE, nrow=dim(mse.nonzero)[1], ncol=dim(mse.nonzero)[2])
for (i in 1:(dim(mse.nonzero)[1])) {msebold.nonzero[i,order(mse.nonzero[i,])[1]] = TRUE}
mseital.nonzero = matrix(FALSE, nrow=dim(mse.nonzero)[1], ncol=dim(mse.nonzero)[2])
for (i in 1:(dim(mse.nonzero)[1])) {mseital.nonzero[i,order(mse.nonzero[i,])[2]] = TRUE}
xtable.printbold(xtable(mse.nonzero, digits=3, align=c('c','c','c','c','c','c')), which.bold=msebold.nonzero, which.ital=mseital.nonzero, include.rownames=FALSE, hline.after=c(0))





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
}
msey2 = cbind(msey, msey.unshrunk, msey.oracular)
colnames(msey2) = c("GWL", "GWL-U", "Oracle")
msey2bold = matrix(FALSE, nrow=dim(msey2)[1], ncol=dim(msey2)[2])
for (i in 1:(dim(msey2)[1])) {msey2bold[i,order(msey2[i,])[1]] = TRUE}
msey2ital = matrix(FALSE, nrow=dim(msey2)[1], ncol=dim(msey2)[2])
for (i in 1:(dim(msey2)[1])) {msey2ital[i,order(msey2[i,])[2]] = TRUE}
xtable.printbold(xtable(msey2, digits=3, caption="Mean squared error of estimates for the response variable $y$ (\\textbf{minimum}, \\emph{next best}).\\label{MSEY}"), which.bold=msey2bold, which.ital=msey2ital, include.rownames=FALSE)

msey = cbind(msey, msey.unshrunk, msey.precon, msey.unshrunk.precon, msey.oracular)
colnames(msey) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")

mseybold = matrix(FALSE, nrow=dim(msey)[1], ncol=dim(msey)[2])
for (i in 1:(dim(msey)[1])) {mseybold[i,order(msey[i,])[1]] = TRUE}

mseyital = matrix(FALSE, nrow=dim(msey)[1], ncol=dim(msey)[2])
for (i in 1:(dim(msey)[1])) {mseyital[i,order(msey[i,])[2]] = TRUE}

xtable.printbold(xtable(msey, digits=3), which.bold=mseybold, which.ital=mseyital, include.rownames=FALSE)




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
}
b22 = cbind(b2, b2.unshrunk, b2.oracular)
colnames(b22) = c("GWL", "GWL-U", "Oracle")
b22bold = matrix(FALSE, nrow=dim(b22)[1], ncol=dim(b22)[2])
for (i in 1:(dim(b22)[1])) {b22bold[i,order(b22[i,])[1]] = TRUE}
b22ital = matrix(FALSE, nrow=dim(b22)[1], ncol=dim(b22)[2])
for (i in 1:(dim(b22)[1])) {b22ital[i,order(b22[i,])[2]] = TRUE}
xtable.printbold(xtable(b22, digits=4, caption="Squared bias of estimates for $\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{bias}"), which.bold=b22bold, which.ital=b22ital, include.rownames=FALSE)

b22.nonzero = cbind(b2.nonzero, b2.unshrunk.nonzero, b2.oracular.nonzero)
colnames(b22.nonzero) = c("GWL", "GWL-U", "Oracle")
b22bold.nonzero = matrix(FALSE, nrow=dim(b22.nonzero)[1], ncol=dim(b22.nonzero)[2])
for (i in 1:(dim(b22.nonzero)[1])) {b22bold.nonzero[i,order(b22.nonzero[i,])[1]] = TRUE}
b22ital.nonzero = matrix(FALSE, nrow=dim(b22.nonzero)[1], ncol=dim(b22.nonzero)[2])
for (i in 1:(dim(b22.nonzero)[1])) {b22ital.nonzero[i,order(b22.nonzero[i,])[2]] = TRUE}
xtable.printbold(xtable(b22.nonzero, digits=4, caption="Squared bias of estimates for $\beta_1$ at locations where $\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{bias-nonzero}"), which.bold=b22bold.nonzero, which.ital=b22ital.nonzero, include.rownames=FALSE)


b2 = cbind(b2, b2.unshrunk, b2.precon, b2.unshrunk.precon, b2.oracular)
colnames(b2) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
b2bold = matrix(FALSE, nrow=dim(b2)[1], ncol=dim(b2)[2])
for (i in 1:(dim(b2)[1])) {b2bold[i,order(b2[i,])[1]] = TRUE}
b2ital = matrix(FALSE, nrow=dim(b2)[1], ncol=dim(b2)[2])
for (i in 1:(dim(b2)[1])) {b2ital[i,order(b2[i,])[2]] = TRUE}
xtable.printbold(xtable(b2, digits=4), which.bold=b2bold, which.ital=b2ital, include.rownames=FALSE)

b2.nonzero = cbind(b2.nonzero, b2.unshrunk.nonzero, b2.precon.nonzero, b2.unshrunk.precon.nonzero, b2.oracular.nonzero)
colnames(b2.nonzero) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
b2bold.nonzero = matrix(FALSE, nrow=dim(b2.nonzero)[1], ncol=dim(b2.nonzero)[2])
for (i in 1:(dim(b2.nonzero)[1])) {b2bold.nonzero[i,order(b2.nonzero[i,])[1]] = TRUE}
b2ital.nonzero = matrix(FALSE, nrow=dim(b2.nonzero)[1], ncol=dim(b2.nonzero)[2])
for (i in 1:(dim(b2.nonzero)[1])) {b2ital.nonzero[i,order(b2.nonzero[i,])[2]] = TRUE}
xtable.printbold(xtable(b2.nonzero, digits=4), which.bold=b2bold.nonzero, which.ital=b2ital.nonzero, include.rownames=FALSE, hline.after=c(0))



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
}
by22 = cbind(by2, by2.unshrunk, by2.oracular)
colnames(by22) = c("GWL", "GWL-U", "Oracle")
by22bold = matrix(FALSE, nrow=dim(by22)[1], ncol=dim(by22)[2])
for (i in 1:(dim(by22)[1])) {by22bold[i,order(by22[i,])[1]] = TRUE}
by22ital = matrix(FALSE, nrow=dim(by22)[1], ncol=dim(by22)[2])
for (i in 1:(dim(by22)[1])) {by22ital[i,order(by22[i,])[2]] = TRUE}
xtable.printbold(xtable(by22, digits=4, caption="Squared bias of estimates for the response variable $y$ (\\textbf{minimum}, \\emph{next best}).\\label{biasy}"), which.bold=by22bold, which.ital=by22ital, include.rownames=FALSE)


by2 = cbind(by2, by2.unshrunk, by2.precon, by2.unshrunk.precon, by2.oracular)
colnames(by2) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")

by2bold = matrix(FALSE, nrow=dim(by2)[1], ncol=dim(by2)[2])
for (i in 1:(dim(by2)[1])) {by2bold[i,order(by2[i,])[1]] = TRUE}

by2ital = matrix(FALSE, nrow=dim(by2)[1], ncol=dim(by2)[2])
for (i in 1:(dim(by2)[1])) {by2ital[i,order(by2[i,])[2]] = TRUE}

xtable.printbold(xtable(by2, digits=4), which.bold=by2bold, which.ital=by2ital, include.rownames=FALSE)



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
}
varx2 = cbind(varx, varx.unshrunk, varx.oracular)
colnames(varx2) = c("GWL", "GWL-U", "Oracle")
varx2bold = matrix(FALSE, nrow=dim(varx2)[1], ncol=dim(varx2)[2])
for (i in 1:(dim(varx2)[1])) {varx2bold[i,order(varx2[i,])[1]] = TRUE}
varx2ital = matrix(FALSE, nrow=dim(varx2)[1], ncol=dim(varx2)[2])
for (i in 1:(dim(varx2)[1])) {varx2ital[i,order(varx2[i,])[2]] = TRUE}
xtable.printbold(xtable(varx2, digits=4, caption="Variance of estimates for $\beta_1$ (\\textbf{minimum}, \\emph{next best}).\\label{varx}"), which.bold=varx2bold, which.ital=varx2ital, include.rownames=FALSE)

varx2.nonzero = cbind(varx.nonzero, varx.unshrunk.nonzero, varx.oracular.nonzero)
colnames(varx2.nonzero) = c("GWL", "GWL-U", "Oracle")
varx2bold.nonzero = matrix(FALSE, nrow=dim(varx2.nonzero)[1], ncol=dim(varx2.nonzero)[2])
for (i in 1:(dim(varx2.nonzero)[1])) {varx2bold.nonzero[i,order(varx2.nonzero[i,])[1]] = TRUE}
varx2ital.nonzero = matrix(FALSE, nrow=dim(varx2.nonzero)[1], ncol=dim(varx2.nonzero)[2])
for (i in 1:(dim(varx2.nonzero)[1])) {varx2ital.nonzero[i,order(varx2.nonzero[i,])[2]] = TRUE}
xtable.printbold(xtable(varx2.nonzero, digits=4, caption="Variance of estimates for $\beta_1$ at locations where $\beta_1 != 0$ (\\textbf{minimum}, \\emph{next best}).\\label{varx-nonzero}"), which.bold=varx2bold.nonzero, which.ital=varx2ital.nonzero, include.rownames=FALSE)


varx = cbind(varx, varx.unshrunk, varx.precon, varx.unshrunk.precon, varx.oracular)
colnames(varx) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
varxbold = matrix(FALSE, nrow=dim(varx)[1], ncol=dim(varx)[2])
for (i in 1:(dim(varx)[1])) {varxbold[i,order(varx[i,])[1]] = TRUE}
varxital = matrix(FALSE, nrow=dim(varx)[1], ncol=dim(varx)[2])
for (i in 1:(dim(varx)[1])) {varxital[i,order(varx[i,])[2]] = TRUE}
xtable.printbold(xtable(varx, digits=4), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE)

varx.nonzero = cbind(varx.nonzero, varx.unshrunk.nonzero, varx.precon.nonzero, varx.unshrunk.precon.nonzero, varx.oracular.nonzero)
colnames(varx.nonzero) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
varxbold.nonzero = matrix(FALSE, nrow=dim(varx.nonzero)[1], ncol=dim(varx.nonzero)[2])
for (i in 1:(dim(varx.nonzero)[1])) {varxbold.nonzero[i,order(varx.nonzero[i,])[1]] = TRUE}
varxital.nonzero = matrix(FALSE, nrow=dim(varx.nonzero)[1], ncol=dim(varx.nonzero)[2])
for (i in 1:(dim(varx.nonzero)[1])) {varxital.nonzero[i,order(varx.nonzero[i,])[2]] = TRUE}
xtable.printbold(xtable(varx.nonzero, digits=4), which.bold=varxbold.nonzero, which.ital=varxital.nonzero, include.rownames=FALSE, hline.after=c(0))



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
}
vary2 = cbind(vary, vary.unshrunk, vary.oracular)
colnames(vary2) = c("GWL", "GWL-U", "Oracle")
vary2bold = matrix(FALSE, nrow=dim(vary2)[1], ncol=dim(vary2)[2])
for (i in 1:(dim(vary2)[1])) {vary2bold[i,order(vary2[i,])[1]] = TRUE}
vary2ital = matrix(FALSE, nrow=dim(vary2)[1], ncol=dim(vary2)[2])
for (i in 1:(dim(vary2)[1])) {vary2ital[i,order(vary2[i,])[2]] = TRUE}
xtable.printbold(xtable(vary2, digits=4, caption="Variance of estimates for the response variable $y$ (\\textbf{minimum}, \\emph{next best}).\\label{vary}"), which.bold=vary2bold, which.ital=vary2ital, include.rownames=FALSE)

vary = cbind(vary, vary.unshrunk, vary.precon, vary.unshrunk.precon, vary.oracular)
colnames(vary) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")
varybold = matrix(FALSE, nrow=dim(vary)[1], ncol=dim(vary)[2])
for (i in 1:(dim(vary)[1])) {varybold[i,order(vary[i,])[1]] = TRUE}
varyital = matrix(FALSE, nrow=dim(vary)[1], ncol=dim(vary)[2])
for (i in 1:(dim(vary)[1])) {varyital[i,order(vary[i,])[2]] = TRUE}
xtable.printbold(xtable(vary, digits=4), which.bold=varybold, which.ital=varyital, include.rownames=FALSE)





errs = list()
for (s in 1:18) {
	err = apply( selection[[s]][['X1']], 2, function(x) {abs((B[['X1']]!=0) - x)})
	vv = c('X2', 'X3', 'X4', 'X5')
	for (v in vv) {
		err = err + selection[[s]][[v]]
	}
	errs[[s]] = err
}

errs.precon = list()
for (s in 1:18) {
	err = apply( selection.precon[[s]][['X1']], 2, function(x) {abs((B[['X1']]!=0) - x)})
	vv = c('X2', 'X3', 'X4', 'X5')
	for (v in vv) {
		err = err + selection.precon[[s]][[v]]
	}
	errs.precon[[s]] = err
}

perfect.selection = data.frame()
for (s in 1:18) {
	perfect.selection = rbind(perfect.selection, c(sum(errs[[s]]==0)/90000, sum(errs.precon[[s]]==0)/90000))
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


print.xtable(xtable(data.frame(tau_x=unique(tau), Frequency=selection.table[['t.x']]), digits=c(0,2,2), caption="Frequency of perfect selection for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0))
print.xtable(xtable(data.frame(tau_sigma=unique(sigma.tau), Frequency=selection.table[['t.e']]), digits=c(0,2,2), caption="Frequency of perfect selection for different settings of $\\tau_{\\sigma}$."), include.rownames=FALSE, hline.after=c(0))
print.xtable(xtable(data.frame(rho=unique(rho), Frequency=selection.table[['rho']]), digits=c(0,2,2), caption="Frequency of perfect selection for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0))


msex.table = list()
msex.table[['t.x']] = cbind(unique(tau), rbind(apply(mse2[t.x[[1]],],2,mean), apply(mse2[t.x[[2]],],2,mean)))
msex.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(mse2[t.e[[1]],],2,mean), apply(mse2[t.e[[2]],],2,mean), apply(mse2[t.e[[3]],],2,mean)))
msex.table[['rho']] = cbind(unique(rho), rbind(apply(mse2[rho.2[[1]],],2,mean), apply(mse2[rho.2[[2]],],2,mean), apply(mse2[rho.2[[3]],],2,mean)))
colnames(msex.table[['t.x']]) = c("tau_x", "GWL", "GWL-U", "Oracle")
colnames(msex.table[['t.e']]) = c("tau_sigma", "GWL", "GWL-U", "Oracle")
colnames(msex.table[['rho']]) = c("rho", "GWL", "GWL-U", "Oracle")

bold.msex = list()
bold.msex[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.msex[['t.x']][i,1+which.min(msex.table[['t.x']][i,2:4])] = TRUE}
bold.msex[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['t.e']][i,1+which.min(msex.table[['t.e']][i,2:4])] = TRUE}
bold.msex[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['rho']][i,1+which.min(msex.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(msex.table[['t.x']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.x']])
xtable.printbold(xtable(msex.table[['t.e']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.e']])
xtable.printbold(xtable(msex.table[['rho']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['rho']])

msex.table.nz = list()
msex.table.nz[['t.x']] = cbind(unique(tau), rbind(apply(mse2.nonzero[t.x[[1]],],2,mean), apply(mse2.nonzero[t.x[[2]],],2,mean)))
msex.table.nz[['t.e']] = cbind(unique(sigma.tau), rbind(apply(mse2.nonzero[t.e[[1]],],2,mean), apply(mse2.nonzero[t.e[[2]],],2,mean), apply(mse2.nonzero[t.e[[3]],],2,mean)))
msex.table.nz[['rho']] = cbind(unique(rho), rbind(apply(mse2.nonzero[rho.2[[1]],],2,mean), apply(mse2.nonzero[rho.2[[2]],],2,mean), apply(mse2.nonzero[rho.2[[3]],],2,mean)))
colnames(msex.table.nz[['t.x']]) = c("tau_x", "GWL", "GWL-U", "Oracle")
colnames(msex.table.nz[['t.e']]) = c("tau_sigma", "GWL", "GWL-U", "Oracle")
colnames(msex.table.nz[['rho']]) = c("rho", "GWL", "GWL-U", "Oracle")

bold.msex.nz = list()
bold.msex.nz[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.msex.nz[['t.x']][i,1+which.min(msex.table[['t.x']][i,2:4])] = TRUE}
bold.msex[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['t.e']][i,1+which.min(msex.table[['t.e']][i,2:4])] = TRUE}
bold.msex[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.msex[['rho']][i,1+which.min(msex.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(msex.table[['t.x']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.x']])
xtable.printbold(xtable(msex.table[['t.e']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['t.e']])
xtable.printbold(xtable(msex.table[['rho']], digits=c(0,2,3,3,3), caption="Mean squared error of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.msex[['rho']])



varx.table = list()
varx.table[['t.x']] = cbind(unique(tau), rbind(apply(varx2[t.x[[1]],],2,mean), apply(varx2[t.x[[2]],],2,mean)))
varx.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(varx2[t.e[[1]],],2,mean), apply(varx2[t.e[[2]],],2,mean), apply(varx2[t.e[[3]],],2,mean)))
varx.table[['rho']] = cbind(unique(rho), rbind(apply(varx2[rho.2[[1]],],2,mean), apply(varx2[rho.2[[2]],],2,mean), apply(varx2[rho.2[[3]],],2,mean)))
colnames(varx.table[['t.x']]) = c("tau_x", "GWL", "GWL-U", "Oracle")
colnames(varx.table[['t.e']]) = c("tau_sigma", "GWL", "GWL-U", "Oracle")
colnames(varx.table[['rho']]) = c("rho", "GWL", "GWL-U", "Oracle")

bold.varx = list()
bold.varx[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.varx[['t.x']][i,1+which.min(varx.table[['t.x']][i,2:4])] = TRUE}
bold.varx[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.varx[['t.e']][i,1+which.min(varx.table[['t.e']][i,2:4])] = TRUE}
bold.varx[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.varx[['rho']][i,1+which.min(varx.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(varx.table[['t.x']], digits=c(0,2,3,3,3), caption="Variance of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.varx[['t.x']])
xtable.printbold(xtable(varx.table[['t.e']], digits=c(0,2,3,3,3), caption="Variance of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.varx[['t.e']])
xtable.printbold(xtable(varx.table[['rho']], digits=c(0,2,3,3,3), caption="Variance of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.varx[['rho']])


bx.table = list()
bx.table[['t.x']] = cbind(unique(tau), rbind(apply(b22[t.x[[1]],],2,mean), apply(b22[t.x[[2]],],2,mean)))
bx.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(b22[t.e[[1]],],2,mean), apply(b22[t.e[[2]],],2,mean), apply(b22[t.e[[3]],],2,mean)))
bx.table[['rho']] = cbind(unique(rho), rbind(apply(b22[rho.2[[1]],],2,mean), apply(b22[rho.2[[2]],],2,mean), apply(b22[rho.2[[3]],],2,mean)))
colnames(bx.table[['t.x']]) = c("tau_x", "GWL", "GWL-U", "Oracle")
colnames(bx.table[['t.e']]) = c("tau_sigma", "GWL", "GWL-U", "Oracle")
colnames(bx.table[['rho']]) = c("rho", "GWL", "GWL-U", "Oracle")

bold.bx = list()
bold.bx[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.bx[['t.x']][i,1+which.min(bx.table[['t.x']][i,2:4])] = TRUE}
bold.bx[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.bx[['t.e']][i,1+which.min(bx.table[['t.e']][i,2:4])] = TRUE}
bold.bx[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.bx[['rho']][i,1+which.min(bx.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(bx.table[['t.x']], digits=c(0,2,3,3,3), caption="Squared bias of the estimate of $\\beta_1$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.bx[['t.x']])
xtable.printbold(xtable(bx.table[['t.e']], digits=c(0,2,3,3,3), caption="Squared bias of the estimate of $\\beta_1$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.bx[['t.e']])
xtable.printbold(xtable(bx.table[['rho']], digits=c(0,2,3,3,3), caption="Squared bias of the estimate of $\\beta_1$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.bx[['rho']])




by.table = list()
by.table[['t.x']] = cbind(unique(tau), rbind(apply(by22[t.x[[1]],],2,mean), apply(by22[t.x[[2]],],2,mean)))
by.table[['t.e']] = cbind(unique(sigma.tau), rbind(apply(by22[t.e[[1]],],2,mean), apply(by22[t.e[[2]],],2,mean), apply(by22[t.e[[3]],],2,mean)))
by.table[['rho']] = cbind(unique(rho), rbind(apply(by22[rho.2[[1]],],2,mean), apply(by22[rho.2[[2]],],2,mean), apply(by22[rho.2[[3]],],2,mean)))
colnames(by.table[['t.x']]) = c("tau_x", "GWL", "GWL-U", "Oracle")
colnames(by.table[['t.e']]) = c("tau_sigma", "GWL", "GWL-U", "Oracle")
colnames(by.table[['rho']]) = c("rho", "GWL", "GWL-U", "Oracle")

bold.by = list()
bold.by[['t.x']] = matrix(FALSE, 2, 4)
for (i in 1:2) {bold.by[['t.x']][i,1+which.min(by.table[['t.x']][i,2:4])] = TRUE}
bold.by[['t.e']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.by[['t.e']][i,1+which.min(by.table[['t.e']][i,2:4])] = TRUE}
bold.by[['rho']] = matrix(FALSE, 3, 4)
for (i in 1:3) {bold.by[['rho']][i,1+which.min(by.table[['rho']][i,2:4])] = TRUE}

xtable.printbold(xtable(by.table[['t.x']], digits=c(0,2,3,3,3), caption="Squared bias of the estimated response variable $y$ for different settings of $\\tau_x$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.by[['t.x']])
xtable.printbold(xtable(by.table[['t.e']], digits=c(0,2,3,3,3), caption="Squared bias of the estimated response variable $y$ for different settings of $\\tau_\\sigma$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.by[['t.e']])
xtable.printbold(xtable(by.table[['rho']], digits=c(0,2,3,3,3), caption="Squared bias of the estimated response variable $y$ for different settings of $\\rho$."), include.rownames=FALSE, hline.after=c(0), which.bold=bold.by[['rho']])



