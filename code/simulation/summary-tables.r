


mse = vector()
mse.oracular = vector()
mse.precon = vector()
mse.unshrunk = vector()
mse.unshrunk.precon = vector()
for (s in 1:18) {
	mse = c(mse, mean(sapply(X1.err[[s]], function(x) {mean(x**2)})))
	mse.oracular = c(mse.oracular, mean(sapply(X1.err.oracular[[s]], function(x) {mean(x**2)})))
	mse.precon = c(mse.precon, mean(sapply(X1.err.precon[[s]], function(x) {mean(x**2)})))
	mse.unshrunk = c(mse.unshrunk, mean(sapply(X1.err.unshrunk[[s]], function(x) {mean(x**2)})))
	mse.unshrunk.precon = c(mse.unshrunk.precon, mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {mean(x**2)})))
}
mse = cbind(mse, mse.unshrunk, mse.precon, mse.unshrunk.precon, mse.oracular)
colnames(mse) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")

msebold = matrix(FALSE, nrow=dim(mse)[1], ncol=dim(mse)[2])
for (i in 1:(dim(mse)[1])) {msebold[i,order(mse[i,])[1]] = TRUE}

mseital = matrix(FALSE, nrow=dim(mse)[1], ncol=dim(mse)[2])
for (i in 1:(dim(mse)[1])) {mseital[i,order(mse[i,])[2]] = TRUE}

xtable.printbold(xtable(mse, digits=3, align=c('c','c','c','c','c','c')), which.bold=msebold, which.ital=mseital, include.rownames=FALSE)


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
for (s in 1:18) {
	b2 = c(b2, mean(sapply(Map(function(z) {mean(sapply(X1.err[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.oracular = c(b2.oracular, mean(sapply(Map(function(z) {mean(sapply(X1.err.oracular[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.precon = c(b2.precon, mean(sapply(Map(function(z) {mean(sapply(X1.err.precon[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.unshrunk = c(b2.unshrunk, mean(sapply(Map(function(z) {mean(sapply(X1.err.unshrunk[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
	b2.unshrunk.precon = c(b2.unshrunk.precon, mean(sapply(Map(function(z) {mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[z]}))**2}, 1:900), identity)))
}
b2 = cbind(b2, b2.unshrunk, b2.precon, b2.unshrunk.precon, b2.oracular)
colnames(b2) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")

b2bold = matrix(FALSE, nrow=dim(b2)[1], ncol=dim(b2)[2])
for (i in 1:(dim(b2)[1])) {b2bold[i,order(b2[i,])[1]] = TRUE}

b2ital = matrix(FALSE, nrow=dim(b2)[1], ncol=dim(b2)[2])
for (i in 1:(dim(b2)[1])) {b2ital[i,order(b2[i,])[1]] = TRUE}

xtable.printbold(xtable(b2, digits=4), which.bold=b2bold, which.ital=b2ital, include.rownames=FALSE)



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
for (s in 1:18) {
	varx = c(varx, mean(sapply(Map(function(z) {var(sapply(X1.err[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.oracular = c(varx.oracular, mean(sapply(Map(function(z) {var(sapply(X1.err.oracular[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.precon = c(varx.precon, mean(sapply(Map(function(z) {var(sapply(X1.err.precon[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.unshrunk = c(varx.unshrunk, mean(sapply(Map(function(z) {var(sapply(X1.err.unshrunk[[s]], function(x) {x[z]}))}, 1:900), identity)))
	varx.unshrunk.precon = c(varx.unshrunk.precon, mean(sapply(Map(function(z) {var(sapply(X1.err.unshrunk.precon[[s]], function(x) {x[z]}))}, 1:900), identity)))
}
varx = cbind(varx, varx.unshrunk, varx.precon, varx.unshrunk.precon, varx.oracular)
colnames(varx) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")

varxbold = matrix(FALSE, nrow=dim(varx)[1], ncol=dim(varx)[2])
for (i in 1:(dim(varx)[1])) {varxbold[i,order(varx[i,])[1]] = TRUE}

varxital = matrix(FALSE, nrow=dim(varx)[1], ncol=dim(varx)[2])
for (i in 1:(dim(varx)[1])) {varxital[i,order(varx[i,])[2]] = TRUE}

xtable.printbold(xtable(varx, digits=4), which.bold=varxbold, which.ital=varxital, include.rownames=FALSE)



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
vary = cbind(vary, vary.unshrunk, vary.precon, vary.unshrunk.precon, vary.oracular)
colnames(by2) = c("AL", "AL-Unshrunk", "AL-Precon", "AL-Precon-Unshrunk", "Oracle")

varybold = matrix(FALSE, nrow=dim(vary)[1], ncol=dim(vary)[2])
for (i in 1:(dim(vary)[1])) {varybold[i,order(vary[i,])[1]] = TRUE}

varyital = matrix(FALSE, nrow=dim(vary)[1], ncol=dim(vary)[2])
for (i in 1:(dim(vary)[1])) {varyital[i,order(vary[i,])[2]] = TRUE}

xtable.printbold(xtable(vary, digits=4), which.bold=varybold, which.ital=vary.ital, include.rownames=FALSE)





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




