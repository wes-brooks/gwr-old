gwlars.ssr = function(bw, x, y, N=N, coords, loc, dist, s, mode.select, verbose, prior.weights, gweight, target, mode, adapt, precondition) {
    gwlars.object = gwlars.fit.inner(x=x, y=y, coords=coords, loc=loc, bw=bw, dist=dist, N=N, s=s, mode.select=mode.select, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode=mode, precondition=precondition)
    resid = sum(gwlars.object[['resid']]**2)
    if (verbose) {cat(paste('loc:(', paste(loc, collapse=","), '), target: ', target, ', bw:', bw, ', ssr:', resid, ', miss:', abs(resid-target), '\n', sep=""))}
    return(abs(resid-target))
}
