mse = vector()
for (i in 1:900) { mse = c(mse, mean(sapply(Y.err[[1]], function(x) {x[i]**2}))) }
mm = matrix(mse, nrow=30, ncol=30)




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

xtable.printbold(xtable(mse, digits=3), which=msebold)


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

xtable.printbold(xtable(msey, digits=3), which=mseybold)



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

xtable.printbold(xtable(perfect.selection, digits=3), which=selectbold)



b = vector()
b.oracular = vector()
b.precon = vector()
b.unshrunk = vector()
b.unshrunk.precon = vector()
for (s in 1:18) {
	b = c(b, mean(sapply(X1.err[[s]], function(x) {mean(x**2)})))
	b.oracular = c(b.oracular, mean(sapply(X1.err.oracular[[s]], function(x) {mean(x**2)})))
	b.precon = c(b.precon, mean(sapply(X1.err.precon[[s]], function(x) {mean(x**2)})))
	b.unshrunk = c(b.unshrunk, mean(sapply(X1.err.unshrunk[[s]], function(x) {mean(x**2)})))
	b.unshrunk.precon = c(b.unshrunk.precon, mean(sapply(X1.err.unshrunk.precon[[s]], function(x) {mean(x**2)})))
}



for s in 1:18 {
	b = vector()
	for (i in 1:900) { b = c(b, mean(sapply(X1.err[[s]], function(x) {x[i]}))) }
	bb = matrix(b, nrow=30, ncol=30)
}

