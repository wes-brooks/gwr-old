gwlars.ssr = function(bw, x, y, coords, loc, dist, s, verbose, prior.weights, gweight, target, mode, adapt) {
    gwlars.object = gwlars.fit.inner(x=x, y=y, coords=coords, loc=loc, bw=bw, dist=dist, s=s, verbose=verbose, gwr.weights=NULL, prior.weights=prior.weights, gweight=gweight, adapt=adapt, mode=mode)
    resid = sum(gwlars.object[['resid']]**2)
    if (verbose) {cat(paste('loc:(', paste(loc, collapse=","), '), target: ', target, ', bw:', bw, ', ssr:', resid, ', miss:', abs(resid-target), '\n', sep=""))}
    return(abs(resid-target))
}
