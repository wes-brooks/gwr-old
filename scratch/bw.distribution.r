bw.distribution = function(trace) {
    spl = smooth.spline(trace, penalty=0)
    xx = range(spl[['x']])
    xxx = seq(xx[1], xx[2], length.out=10001)
    smoothed = predict(spl, xxx)
    
    cbind(xxx, cumsum(exp(-smoothed[['y']]))/sum(exp(-smoothed[['y']])))
}




n = nrow(trace)
d1 = diag(rep(-1,n))
d1[2:n,1:(n-1)] = d1[2:n,1:(n-1)] + diag(rep(1,n-1))
diff1 = t(Matrix(trace)) %*% d1
diff1 = diff1[,1:(n-1)]


d2 = diag(rep(-1,n-1))
d2[2:(n-1),1:(n-2)] = d2[2:(n-1),1:(n-2)] + diag(rep(1,n-2))
diff2 = diff1 %*% d2
diff2 = t(diff2[,1:(n-2)])
