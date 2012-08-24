\name{gwglmnet.sel}
\alias{gwglmnet.sel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bandwidth selection in a GW-GLM model (bandwidth in terms of nearest neighbors or distance).
}
\description{
Bandwidth selection in a GW-GLM model (bandwidth in terms of nearest neighbors or distance).
}
\usage{
gwglmnet.sel(formula, data = list(), coords, adapt = FALSE, nearest.neighbors = FALSE, gweight = gwr.Gauss, s, method = "cv", verbose = FALSE, longlat = FALSE, family, weights, tol = .Machine$double.eps^0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{adapt}{
%%     ~~Describe \code{adapt} here~~
}
  \item{nearest.neighbors}{
%%     ~~Describe \code{nearest.neighbors} here~~
}
  \item{gweight}{
%%     ~~Describe \code{gweight} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{longlat}{
%%     ~~Describe \code{longlat} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wesley Brooks
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, data = list(), coords, adapt = FALSE, nearest.neighbors = FALSE, 
    gweight = gwr.Gauss, s, method = "cv", verbose = FALSE, longlat = FALSE, 
    family, weights, tol = .Machine$double.eps^0.25) 
{
    if (!is.logical(adapt)) 
        stop("adapt must be logical")
    if (is(data, "Spatial")) {
        if (!missing(coords)) 
            warning("data is Spatial* object, ignoring coords argument")
        coords <- coordinates(data)
        if ((is.null(longlat) || !is.logical(longlat)) && !is.na(is.projected(data)) && 
            !is.projected(data)) {
            longlat <- TRUE
        }
        else longlat <- FALSE
        data <- as(data, "data.frame")
    }
    if (is.null(longlat) || !is.logical(longlat)) 
        longlat <- FALSE
    if (missing(coords)) 
        stop("Observation coordinates have to be given")
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
    weights <- as.vector(model.extract(mf, "weights"))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) 
        weights <- rep(as.numeric(1), dp.n)
    if (any(is.na(weights))) 
        stop("NAs in weights")
    if (any(weights < 0)) 
        stop("negative weights")
    y <- model.extract(mf, "response")
    x <- model.matrix(mt, mf)
    n = dim(coords)[1]
    if (longlat) {
        D = as.matrix(earth.dist(coords), n, n)
    }
    else {
        Xmat = matrix(rep(coords[, 1], times = n), n, n)
        Ymat = matrix(rep(coords[, 2], times = n), n, n)
        D = sqrt((Xmat - t(Xmat))^2 + (Ymat - t(Ymat))^2)
    }
    if (nearest.neighbors) {
        beta1 <- 0
        beta2 <- 1
    }
    else {
        bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
        difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
        if (any(!is.finite(difmin))) 
            difmin[which(!is.finite(difmin))] <- 0
        beta1 <- difmin/1000
        beta2 <- 2 * difmin
    }
    opt <- optimize(gwglmnet.cv.f, lower = beta1, upper = beta2, 
        maximum = FALSE, tol = tol, formula = formula, coords = coords, 
        s = s, gweight = gweight, verbose = verbose, longlat = longlat, 
        data = data, D = D, weights = weights, adapt = adapt, 
        nn = nearest.neighbors, family = family)
    bdwt <- opt$minimum
    res <- bdwt
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
