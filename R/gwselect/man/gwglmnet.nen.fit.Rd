\name{gwglmnet.nen.fit}
\alias{gwglmnet.nen.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a GW-GLM model using Nearest Effective Neighbors for bandwidth selection.
}
\description{
Fit a GW-GLM model using Nearest Effective Neighbors for bandwidth selection.
}
\usage{
gwglmnet.nen.fit(x, y, coords, D, s = NULL, verbose, family, prior.weights, gweight, bw, beta1, beta2, type = "pearson", tol = 1e-25, longlat = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{coords}{
%%     ~~Describe \code{coords} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{prior.weights}{
%%     ~~Describe \code{prior.weights} here~~
}
  \item{gweight}{
%%     ~~Describe \code{gweight} here~~
}
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{beta1}{
%%     ~~Describe \code{beta1} here~~
}
  \item{beta2}{
%%     ~~Describe \code{beta2} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
}
  \item{longlat}{
%%     ~~Describe \code{longlat} here~~
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
function (x, y, coords, D, s = NULL, verbose, family, prior.weights, 
    gweight, bw, beta1, beta2, type = "pearson", tol = 1e-25, 
    longlat = FALSE) 
{
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    gwglmnet.object = list()
    cv.error = list()
    for (i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[, 1] == coords.unique[i, 1] & 
            coords[, 2] == coords.unique[i, 2])
        dist = D[i, ]
        bandwidth = optimize(gwglmnet.ssr, lower = beta1, upper = beta2, 
            maximum = FALSE, tol = bw/10, x = x, y = y, colocated = colocated, 
            s = s, gweight = gweight, verbose = verbose, dist = dist, 
            prior.weights = prior.weights, family = family, target = bw, 
            type = type)$minimum
        cat(paste("For i=", i, ", bw=", bandwidth, ".\n", sep = ""))
        weight.matrix = gweight(D, bandwidth)
        loow = weight.matrix[i, -colocated]
        prior.loow = prior.weights[-colocated]
        reps = length(colocated)
        w <- prior.loow * loow
        if (sum(loow) == 0) {
            return(list(cv.error = Inf))
        }
        reps = length(colocated)
        xx = as.matrix(x[-colocated, ])
        yy = as.matrix(y[-colocated])
        if (family == "binomial" && (abs(sum(yy * w) - sum(w)) < 
            1e-04 || sum(yy * w) < 1e-04)) {
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yy * 
                w), ", sum of weights=", sum(w), "\n", sep = ""))
            model[[i]] = NULL
            cv.error[[i]] = 0
            s.optimal = c(s.optimal, max(s))
        }
        else if (family == "binomial") {
            model[[i]] = glmnet(x = xx, y = cbind(1 - yy, yy), 
                weights = w, family = family, lambda = s)
            predictions = predict(model[[i]], newx = matrix(x[colocated, 
                ], nrow = reps, ncol = dim(xx)[2]), s = s, type = "response")
            cv.error[[i]] = colSums(abs(matrix(predictions - 
                matrix(y[colocated], nrow = reps, ncol = length(s)), 
                nrow = reps, ncol = length(s))))
            s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
        }
        else {
            model[[i]] = glmnet(x = xx, y = yy, weights = w, 
                family = family, lambda = s)
            predictions = predict(model[[i]], newx = matrix(x[colocated, 
                ], nrow = reps, ncol = dim(xx)[2]), s = s, type = "response")
            cv.error[[i]] = colSums(abs(matrix(predictions - 
                matrix(y[colocated], nrow = reps, ncol = length(s)), 
                nrow = reps, ncol = length(s))))
            s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
        }
        if (verbose) {
            cat(paste(i, "\n", sep = ""))
        }
    }
    gwglmnet.object[["coef.scale"]] = NULL
    gwglmnet.object[["model"]] = model
    gwglmnet.object[["s"]] = s.optimal
    gwglmnet.object[["mode"]] = mode
    gwglmnet.object[["coords"]] = coords.unique
    gwglmnet.object[["cv.error"]] = cv.error
    gwglmnet.object[["s.range"]] = s
    class(gwglmnet.object) = "gwglmnet.object"
    return(gwglmnet.object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
