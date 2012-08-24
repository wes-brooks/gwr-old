neighbor.diff = function(bw, D, q.target, weight.function, verbose) {
    if (verbose) {cat(paste("bandwidth: ", bw, '\n', sep=''))}

    nobs = length(D)
    W = weight.function(D, bw)
    q.this = sum(W)/nobs

    if (verbose) {cat(paste("total weight: ", q.this, '\n', sep=''))}

    return(abs(q.target - q.this))
}
