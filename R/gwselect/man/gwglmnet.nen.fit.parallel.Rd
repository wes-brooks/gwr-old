\name{gwglmnet.nen.fit.parallel}
\alias{gwglmnet.nen.fit.parallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multicore-aware function to fit a GW-GLM model using the LASSO for variable selection and Nearest Effective Neighbors for bandwidth selection.
}
\description{
Multicore-aware function to fit a GW-GLM model using the LASSO for variable selection and Nearest Effective Neighbors for bandwidth selection.
}
\usage{
gwglmnet.nen.fit.parallel(x, y, coords, D, s, verbose, family, prior.weights, gweight, target, beta1, beta2, type = "pearson", tol = 1e-25, longlat = FALSE)
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
  \item{target}{
%%     ~~Describe \code{target} here~~
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
function (x, y, coords, D, s, verbose, family, prior.weights, 
    gweight, target, beta1, beta2, type = "pearson", tol = 1e-25, 
    longlat = FALSE) 
{
    coords.unique = unique(coords)
    n = dim(coords.unique)[1]
    gwglmnet.object = foreach(i = 1:n, .packages = "glmnet", 
        .errorhandling = "remove") \%dopar\% {
        colocated = which(coords[, 1] == coords.unique[i, 1] & 
            coords[, 2] == coords.unique[i, 2])
        dist = D[i, ]
        opt = optimize(gwglmnet.ssr, lower = beta1, upper = beta2, 
            maximum = FALSE, tol = target/1000, x = x, y = y, 
            colocated = colocated, s = s, gweight = gweight, 
            verbose = verbose, dist = dist, prior.weights = prior.weights, 
            family = family, target = target, type = type)
        bandwidth = opt$minimum
        cat(paste("For i=", i, ", bw=", bandwidth, ", tolerance=", 
            target/1000, ", miss=", sqrt(opt$objective), ".\n", 
            sep = ""))
        loow = gweight(D[i, -colocated], bandwidth)
        prior.loow = prior.weights[-colocated]
        w <- prior.loow * loow
        reps = length(colocated)
        if (sum(loow) == 0) {
            return(list(cv.error = Inf))
        }
        xx = as.matrix(x[-colocated, ])
        yy = as.matrix(y[-colocated])
        if (family == "binomial" && (abs(sum(yy * w) - sum(w)) < 
            1e-04 || sum(yy * w) < 1e-04)) {
            cat(paste("Abort. i=", i, ", weighted sum=", sum(yy * 
                w), ", sum of weights=", sum(w), "\n", sep = ""))
            model = NULL
            cv.error = 0
            s.optimal = max(s)
        }
        else if (family == "binomial") {
            model = glmnet(x = xx, y = cbind(1 - yy, yy), weights = w, 
                family = family, lambda = s)
            predictions = predict(model, newx = matrix(x[colocated, 
                ], nrow = reps, ncol = dim(xx)[2]), s = s, type = "response")
            cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], 
                nrow = reps, ncol = length(s)), nrow = reps, 
                ncol = length(s))))
            s.optimal = s[which.min(cv.error)]
        }
        else {
            model = glmnet(x = xx, y = yy, weights = w, family = family, 
                lambda = s)
            ll = model$lambda
            predictions = predict(model, newx = matrix(x[colocated, 
                ], nrow = reps, ncol = dim(xx)[2]), s = ll, type = "response")
            cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], 
                nrow = reps, ncol = length(ll)), nrow = reps, 
                ncol = length(ll))))
            s.optimal = ll[which.min(cv.error)]
        }
        if (verbose) {
            cat(paste(i, "\n", sep = ""))
        }
        list(model = model, cv.error = cv.error, s = s.optimal, 
            index = i)
    }
    print("returning from gwglmnet.nen.fit.parallel")
    class(gwglmnet.object) = "gwglmnet.object"
    return(gwglmnet.object)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
