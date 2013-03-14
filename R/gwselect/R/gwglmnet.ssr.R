gwglmnet.ssr = function(bw, x, y, family, coords, loc, dist, s, verbose, prior.weights, gweight, target, adapt, precondition, interact, alpha) {
    gwglmnet.object = gwglmnet.fit.inner(x=x, y=y, family=family, coords=coords, loc=loc, bw=bw, dist=dist, s=s, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, precondition=precondition, interact=interact, alpha=alpha)
    loss = gwglmnet.object[['loss']]
    if (verbose) {cat(paste('loc:(', paste(loc, collapse=","), '), target: ', target, ', bw:', bw, ', ssr:', loss, ', miss:', abs(loss-target), '\n', sep=""))}
    return(abs(loss-target))
}
