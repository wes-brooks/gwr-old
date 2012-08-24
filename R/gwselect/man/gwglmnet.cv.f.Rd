\name{gwglmnet.cv.f}
\alias{gwglmnet.cv.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform cross-validation for bandwidth selection in a gw-glm model.
}
\description{
Perform cross-validation for bandwidth selection in a gw-glm model.
}
\usage{
gwglmnet.cv.f(formula, data, bw, coords, gweight, verbose, adapt, longlat, s, family, weights, nn, D = NULL, ...)
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
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{nn}{
%%     ~~Describe \code{nn} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
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
    longlat, s, family, weights, nn, D = NULL, ...) 
{
    cat(paste("Beginning with bandwidth: ", bw, "\n", sep = ""))
    gwglmnet.model = gwglmnet(formula = formula, data = data, 
        coords = coords, gweight = gweight, bw = bw, verbose = verbose, 
        longlat = longlat, adapt = adapt, s = s, family = family, 
        weights = weights, nearest.neighbors = nn, D = D)
    cv.error = sum(sapply(gwglmnet.model[["cv.error"]], min))
    cat(paste("Bandwidth: ", bw, ". CV error: ", cv.error, "\n", 
        sep = ""))
    return(cv.error)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
