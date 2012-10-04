\name{registerCores}
\alias{registerCores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Register multiple cores for parallelization via doMC
}
\description{
Register multiple cores for parallelization via doMC
}
\usage{
gwglmnet(formula, data, coords, gweight, bw, D = NULL, verbose = FALSE, longlat = FALSE, adapt = FALSE, s, family, weights = NULL, nearest.neighbors = FALSE)
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
  \item{gweight}{
%%     ~~Describe \code{gweight} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{longlat}{
%%     ~~Describe \code{longlat} here~~
}
  \item{adapt}{
%%     ~~Describe \code{adapt} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{nearest.neighbors}{
%%     ~~Describe \code{nearest.neighbors} here~~
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
function (formula, data, coords, gweight, bw, D = NULL, verbose = FALSE, 
    longlat = FALSE, adapt = FALSE, s, family, weights = NULL, 
    nearest.neighbors = FALSE) 
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
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    dp.n <- length(model.extract(mf, "response"))
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
    if (is.null(D)) {
        n = dim(coords)[1]
        if (longlat) {
            D = as.matrix(earth.dist(coords), n, n)
        }
        else {
            Xmat = matrix(rep(coords[, 1], times = n), n, n)
            Ymat = matrix(rep(coords[, 2], times = n), n, n)
            D = sqrt((Xmat - t(Xmat))^2 + (Ymat - t(Ymat))^2)
        }
    }
    if (!nearest.neighbors) {
        weight.matrix = gweight(D, bw)
    }
    else {
        n = dim(D)[1]
        bandwidths = sapply(1:n, function(x) {
            neighbor.weight(q = bw, D = D[x, ], weight.function = gweight, 
                verbose = verbose, tol = 0.001)
        })
        weight.matrix = as.matrix(rbind(sapply(1:n, function(k) {
            gweight(as.vector(D[k, ]), as.numeric(bandwidths[1, 
                k]))
        })), n, n)
    }
    if (!adapt) {
        res = gwglmnet.fit(x, y, coords, weight.matrix, s, verbose, 
            family, weights)
    }
    else {
        res = gwglmnet.adaptive.fit(x, y, coords, weight.matrix, 
            s, verbose, family, weights)
    }
    res[["data"]] = data
    res[["response"]] = as.character(formula[[2]])
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
