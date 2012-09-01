gwlars.knn = function(bw, x, y, coords, loc, dist, s, verbose, prior.weights, total.weight, gweight, target, adapt, precondition) {
    gwr.weights = gweight(dist, bw)
    w = gwr.weights * prior.weights
    prop = sum(w)/total.weight

    if (verbose) {cat(paste('loc:(', paste(loc, collapse=","), '), target:', target, ', bw:', bw, ', proportion of weight:', prop, ', miss:', abs(prop-target), '\n', sep=""))}
    return(abs(prop-target))
}
