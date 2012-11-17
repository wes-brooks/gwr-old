vars = c("X1", "X2", "X3", "X4", "X5")
xx = c(0, 30)
yy = c(0,1)

plot(apply(matrix(ss[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="selection frequency")
for (k in 2:length(vars)) {
    v = vars[k]
    par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
    plot(apply(matrix(ss[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
}
legend(vars, lty=1:length(vars), bty='n', x='topleft')


dev.new()
vars = c("X1", "X2", "X3", "X4", "X5", "(Intercept)")
plot(apply(matrix(cs[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
for (k in 2:length(vars)) {
    v = vars[k]
    par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
    plot(apply(matrix(cs[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
}
abline(h=0.95, lty=3, col='red')
legend(vars, lty=1:length(vars), bty='n', x='bottomleft')


dev.new()
plot(apply(matrix(cb[['X1']], 30, 30),1,mean), type='l', lty=1, bty='n', xlim=xx, ylim=yy, xlab="location", ylab="95% CI coverage frequency")
for (k in 2:length(vars)) {
    v = vars[k]
    par(new=TRUE, xaxt='n', yaxt='n', ann=FALSE, bty='n')
    plot(apply(matrix(cb[[v]], 30, 30),1,mean), type='l', lty=k, xlim=xx, ylim=yy)
}
abline(h=0.95, lty=3, col='red')
legend(vars, lty=1:length(vars), bty='n', x='bottomleft')


coverage.b = data.frame()
coverage.s = data.frame()
selection = data.frame()

cbv = vector()
csv = vector()
sv = vector()
for (k in 1:length(vars)) {
    v = vars[k]
    csv = c(csv, mean(cs[[v]]))
    cbv = c(cbv, mean(cb[[v]]))
    
    if (v!="(Intercept)") {
        sv = c(sv, mean(ifelse(B[[v]]==0, 1-ss[[v]], ss[[v]])))
    }
}

coverage.b = rbind(coverage.b, cbv)
coverage.s = rbind(coverage.s, csv)
selection = rbind(selection, sv)

colnames(coverage.b) = vars
colnames(coverage.s) = vars
colnames(selection) = vars[1:5]
