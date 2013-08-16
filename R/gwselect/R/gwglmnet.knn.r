gwglmnet.knn = function(bw, coords, indx=NULL, loc, dist, verbose, prior.weights, total.weight, gweight, target) {
    gwr.weights = gweight(dist, bw)

    if (!is.null(indx)) {
        gwr.weights = gwr.weights[indx]
    }

    w = gwr.weights * prior.weights
    prop = sum(w)/total.weight

    #if (verbose) {cat(paste('loc:(', paste(loc, collapse=","), '), target:', target, ', bw:', bw, ', proportion of weight:', prop, ', miss:', abs(prop-target), '\n', sep=""))}
    return(abs(prop-target))
}
