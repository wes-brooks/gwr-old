\name{gwglmnet.nen.cv.f}
\alias{gwglmnet.nen.cv.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation for selection of tuning parameter in a GW-GLM model using Nearest Effective Neighbors for bandwidth selection.
}
\description{
Cross-validation for selection of tuning parameter in a GW-GLM model using Nearest Effective Neighbors for bandwidth selection.
}
\usage{
gwglmnet.nen.cv.f(formula, data, bw, coords, gweight, verbose, adapt, longlat, s = NULL, beta1, beta2, family, weights = NULL, D = NULL, tolerance = .Machine$double.eps^0.25, type = "pearson", parallel = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{gweight}{
%%     ~~Describe \code{gweight} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{adapt}{
%%     ~~Describe \code{adapt} here~~
}
  \item{longlat}{
%%     ~~Describe \code{longlat} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
}
  \item{beta2}{
%%     ~~Describe \code{beta2} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{tolerance}{
%%     ~~Describe \code{tolerance} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (formula, data, bw, coords, gweight, verbose, adapt, 
    longlat, s = NULL, beta1, beta2, family, weights = NULL, 
    D = NULL, tolerance = .Machine$double.eps^0.25, type = "pearson", 
    parallel = FALSE, ...) 
{
    cat(paste("Beginning with target SSR: ", bw, ", tolerance: ", 
        tolerance, "\n", sep = ""))
    gwglmnet.model = gwglmnet.nen(formula = formula, data = data, 
        coords = coords, gweight = gweight, bw = bw, verbose = verbose, 
        longlat = longlat, adapt = adapt, s = s, family = family, 
        weights = weights, D = D, tol = tolerance, beta1 = beta1, 
        beta2 = beta2, type, parallel = parallel)
    print(gwglmnet.model[["model"]][["cv.error"]])
    print(names(gwglmnet.model))
    print(gwglmnet.model[["model"]])
    cv.error = sum(sapply(gwglmnet.model[["model"]], function(x) min(x[["cv.error"]])))
    cat(paste("Bandwidth: ", bw, ". CV error: ", cv.error, "\n", 
        sep = ""))
    return(cv.error)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
