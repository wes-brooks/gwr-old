\name{gwglmnet.nen.sel}
\alias{gwglmnet.nen.sel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bandwidth selection using Nearest Effective Neighbors in a GW-GLM model.
}
\description{
Bandwidth selection using Nearest Effective Neighbors in a GW-GLM model.
}
\usage{
gwglmnet.nen.sel(formula, data = list(), coords, adapt = FALSE, gweight = gwr.Gauss, s = NULL, method = "cv", verbose = FALSE, longlat = FALSE, family, weights = NULL, tol = .Machine$double.eps^0.25, type, parallel = FALSE)
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
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (formula, data = list(), coords, adapt = FALSE, gweight = gwr.Gauss, 
    s = NULL, method = "cv", verbose = FALSE, longlat = FALSE, 
    family, weights = NULL, tol = .Machine$double.eps^0.25, type, 
    parallel = FALSE) 
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
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) 
        weights <- rep(1, dp.n)
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
    model = glm(formula = formula, data = data, family = family, 
        weights = weights)
    SSR = sum((weights * residuals(model, type = type))^2)
    cat(paste("The SSR from the global model is: ", SSR, "\n", 
        sep = ""))
    nloc = unique(coords)
    lowerSSR <- SSR/5000
    upperSSR <- SSR
    bbox <- cbind(range(coords[, 1]), range(coords[, 2]))
    difmin <- spDistsN1(bbox, bbox[2, ], longlat)[1]
    if (any(!is.finite(difmin))) 
        difmin[which(!is.finite(difmin))] <- 0
    beta1 <- difmin/1000
    beta2 <- difmin
    cat(paste("Maximum distance: ", difmin, "\n", sep = ""))
    opt <- optimize(gwglmnet.nen.cv.f, lower = lowerSSR, upper = upperSSR, 
        maximum = FALSE, tol = tol, tolerance = tol, formula = formula, 
        coords = coords, s = s, beta1 = beta1, beta2 = beta2, 
        gweight = gweight, verbose = verbose, longlat = longlat, 
        data = data, D = D, weights = weights, adapt = adapt, 
        family = family, type = type, parallel = parallel)
    bdwt <- opt$minimum
    res <- bdwt
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
