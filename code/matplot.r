matplot <- function (x, cs1 = c(0, 1), cs2 = c(0, 1), cs3 = c(0, 1), extremes = NA, 
    cellcolors = NA, show.legend = FALSE, nslices = 10, xlab = "Column", 
    ylab = "Row", do.hex = FALSE, axes = TRUE, show.values = FALSE, 
    vcol = NA, vcex = 1, border = "black", na.color = NA, xrange = NULL, 
    color.spec = "rgb", yrev = TRUE, xat = NULL, yat = NULL, 
    Hinton = FALSE, ...) 
{
    if (is.matrix(x) || is.data.frame(x))
    {
        xdim <- dim(x)
        
        if(is.vector(xrange))
        {  
            if(xrange[1] > min(x, na.rm=TRUE) || xrange[2] < max(x, na.rm=TRUE) || length(xrange) != 2)
            {
                xrange = vector()
                xrange[1] = min(x, na.rm=TRUE) - min(x, na.rm=TRUE) %% 0.1                
                xrange[2] = max(x, na.rm=TRUE) + 0.1 - (max(x, na.rm=TRUE) %% 0.1)
            }
        }
        else
        {
            xrange = vector()
            xrange[1] = min(x, na.rm=TRUE) - min(x, na.rm=TRUE) %% 0.1                
            xrange[2] = max(x, na.rm=TRUE) + 0.1 - (max(x, na.rm=TRUE) %% 0.1)   
        }
        if (is.data.frame(x)) 
            x <- unlist(x)
        else x <- as.vector(x)
        oldpar <- par("xaxs", "yaxs", "xpd", "mar")
        par(xaxs="i", yaxs="i")
        if(do.hex) 
            par(mar = c(5, 4, 4, 4))
        plot(c(0, xdim[2]), c(0, xdim[1]), xlab = xlab, ylab = ylab, 
            type = "n", axes = FALSE, ...)
        oldpar$usr <- par("usr")
        if(!do.hex)
        {
            box()
            pos <- 0
        }
        else pos <- -0.3
        if(axes)
        {
            if (is.null(xat)) 
                xat <- pretty(0:xdim[2])[-1]
            axis(1, at = xat - 0.5, labels = xat, pos = pos)
            if (is.null(yat)) 
                yat <- pretty(0:xdim[1])[-1]
            axis(2, at = xdim[1] - yat + 0.5, labels = yat)
        }
        if (all(is.na(cellcolors)))
        {
            if (Hinton)
            {
                if (is.na(extremes[1])) 
                  extremes <- c("black", "white")
                cellcolors <- extremes[(x > 0) + 1]
            }


            else cellcolors <- color.scale(x, cs1, cs2, cs3, extremes=extremes, na.color=na.color, color.spec=color.spec, xrange=xrange)
        }
        if (is.na(vcol)) 
            vcol <- ifelse(colSums(col2rgb(cellcolors) * c(1, 
                1.4, 0.6)) < 350, "white", "black")
        if (Hinton)
        {
            if (any(x < 0 | x > 1)) 
                cellsize <- matrix(rescale(abs(x), c(0, 1)), 
                  nrow = 10)
        }
        else cellsize <- matrix(1, nrow = xdim[1], ncol = xdim[2])
        if (do.hex)
        {
            par(xpd = TRUE)
            offset <- 0
            if (length(border) < xdim[1] * xdim[2]) 
                border <- rep(border, length.out = xdim[1] * 
                  xdim[2])
            for (row in 1:xdim[1])
            {
                for (column in 0:(xdim[2] - 1))
                {
                  hexagon(column+offset, xdim[1]-row, unitcell=cellsize[row, column+1], col=cellcolors[row + xdim[1]*column], border=border[row + xdim[1]*column])
                  if (show.values) 
                    text(column+offset+0.5, xdim[1]-row+0.5, x[row+column*xdim[1]], col=vcol[row + xdim[1]*column], cex=vcex)
                }
                offset <- ifelse(offset, 0, 0.5)
            }
            par(xpd=FALSE)
        }
        else
        {
            if (Hinton) 
                inset <- (1 - cellsize)/2
            else inset <- 0
            if (yrev)
            {
                y0 <- rep(seq(xdim[1] - 1, 0, by = -1), xdim[2]) + inset
                y1 <- rep(seq(xdim[1], 1, by = -1), xdim[2]) - inset
            }
            else
            {
                y0 <- rep(0:(xdim[1] - 1), xdim[2]) + inset
                y1 <- rep(1:xdim[1], xdim[2]) - inset
            }
            rect(sort(rep((1:xdim[2]) - 1, xdim[1])) + inset, 
                y0, sort(rep(1:xdim[2], xdim[1])) - inset, y1, 
                col = cellcolors, border = border)
            if (show.values)
            {
                if (yrev) texty <- rep(seq(xdim[1] - 0.5, 0, by = -1), xdim[2])
                else texty <- rep(seq(0.5, xdim[1] - 0.5, by = 1), xdim[2])
                text(sort(rep((1:xdim[2]) - 0.5, xdim[1])), texty, 
                  round(x, show.values), col = vcol, cex = vcex)
            }
        }
        naxs <- which(is.na(x))
        xy <- par("usr")
        plot.din <- par("din")
        plot.pin <- par("pin")
        bottom.gap <- (xy[3] - xy[4]) * (plot.din[2] - plot.pin[2])/(2 * plot.pin[2])
        grx1 <- xy[1]
        gry1 <- bottom.gap * 0.95
        grx2 <- xy[1] + (xy[2] - xy[1])/4
        gry2 <- bottom.gap * 0.8
        if (length(cellcolors) > 1)
        {
            colmat <- col2rgb(color.scale(xrange, cs1, cs2, cs3, extremes=extremes, na.color=na.color, color.spec=color.spec, xrange=xrange))
            cs1 <- colmat[1, ]/255
            cs2 <- colmat[2, ]/255
            cs3 <- colmat[3, ]/255
            color.spec <- "rgb"
        }
        rect.col <- color.scale(1:nslices, cs1, cs2, cs3, color.spec = color.spec)
        if (show.legend) 
            color.legend(grx1, gry1, grx2, gry2, round(xrange, show.legend), rect.col=rect.col)
        par(oldpar)
    }
    else cat("x must be a data frame or matrix\n")
}
