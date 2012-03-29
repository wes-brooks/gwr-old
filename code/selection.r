r = function(x, xi) { sqrt(sum((x-xi)**2)) }
r2 = function(x, y, xi, yi) { return( sqrt((x-xi)**2 + (y-yi)**2) ) }


eta = function(x, xi)
{
    if (r(x, xi)==0) return(0)
    else return(r(x, xi)**2 * log(r(x, xi)))
}


eta.i = function(x, y, xi, yi)
{
    if(x==xi & y==yi)
        { return(0) }
    else
        { return(r2(x, y, xi, yi)**2 * log(r2(x, y, xi, yi))) }
}


eta.ii = function(x, y, xi, yi)
{
    return( ifelse(x==xi & y==yi, 0, r2(x, y, xi, yi)**2 * log(r2(x, y, xi, yi))) )
}


xlim = c(0,100)
ylim = c(0,100)

xx = seq(xlim[1], xlim[2]) * 4 / (xlim[2] - xlim[1])
yy = seq(ylim[1], ylim[2]) * 4 / (ylim[2] - ylim[1])

mat = matrix(nrow=length(yy), ncol=length(xx))

for (i in 1:length(xx)) {
    for (j in 1:length(yy)) {
        mat[j,i] = eta(c(xx[j], yy[i]), c(2, 2))
    }
}



E = function(y, xi, yi) { integrate(eta.ii, lower=0, upper=10, y=y, xi=xi, yi=yi) }

integrate(E, lower=0, upper=10, xi=5, yi=5)


