heat_leg <- function(range, r=c(0,1), g=0, b=c(1,0), slices=10)
{
    testcol<-color.gradient(r, g, b, nslices=slices)
    color.legend(0,-4,20,-2, "",testcol)
    mtext(range[1], side=1, line=2, outer=FALSE, cex=1, at=0)
    mtext(range[2], side=1, line=2, outer=FALSE, cex=1, at=20)
}