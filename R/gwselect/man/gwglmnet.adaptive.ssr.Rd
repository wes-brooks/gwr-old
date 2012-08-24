\name{gwglmnet.adaptive.ssr}
\alias{gwglmnet.adaptive.ssr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get the sum of squared residuals in for a geographically-weighted GLM
}
\description{
Get the sum of squared residuals in for a geographically-weighted GLM
}
\usage{
gwglmnet.adaptive.ssr(bw, x, y, colocated, dist, s, verbose, family, prior.weights, gweight, type, target)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bw}{
%%     ~~Describe \code{bw} here~~
}
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{colocated}{
%%     ~~Describe \code{colocated} here~~
}
  \item{dist}{
%%     ~~Describe \code{dist} here~~
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
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
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
function (bw, x, y, colocated, dist, s, verbose, family, prior.weights, 
    gweight, type, target) 
{
    reps = length(colocated)
    loow = gweight(dist, bw)[-colocated]
    w <- prior.weights[-colocated] * loow
    xx = as.matrix(x[-colocated, ])
    yy = as.matrix(y[-colocated])
    m <- ncol(xx)
    n <- nrow(xx)
    one <- rep(1, n)
    meanx <- drop(one \%*\% xx)/n
    x.centered <- scale(xx, meanx, FALSE)
    normx <- sqrt(drop(one \%*\% (x.centered^2)))
    names(normx) <- NULL
    xs = x.centered
    for (k in 1:dim(x.centered)[2]) {
        if (normx[k] != 0) {
            xs[, k] = xs[, k]/normx[k]
        }
        else {
            xs[, k] = rep(0, dim(xs)[1])
            normx[k] = Inf
        }
    }
    glm.step = try(glm(yy ~ xs, family = family, weights = w))
    if (class(glm.step) == "try-error") {
        cat(paste("Couldn't make a model for finding the SSR at bandwidth ", 
            bw, "\n", sep = ""))
        return(Inf)
    }
    beta.glm = glm.step$coeff[2:(m + 1)]
    adapt.weight = abs(beta.glm)
    for (k in 1:dim(x.centered)[2]) {
        if (!is.na(adapt.weight[k])) {
            xs[, k] = xs[, k] * adapt.weight[k]
        }
        else {
            xs[, k] = rep(0, dim(xs)[1])
            adapt.weight[k] = 0
        }
    }
    if (family == "binomial") {
        model = glmnet(x = xs, y = cbind(1 - yy, yy), weights = w, 
            family = family, lambda = s)
    }
    else {
        model = glmnet(x = xs, y = yy, weights = w, family = family, 
            lambda = s)
    }
    ll = model$lambda
    xs.colocated = (x[colocated, ] - meanx) * adapt.weight/normx
    predictions = predict(model, newx = matrix(xs.colocated, 
        nrow = reps, ncol = dim(xs)[2]), s = ll, type = "response", 
        )
    cv.error = colSums(abs(matrix(predictions - matrix(y[colocated], 
        nrow = reps, ncol = length(ll)), nrow = reps, ncol = length(ll))))
    s.optimal = ll[which.min(cv.error)]
    fitted = predict(model, newx = xs, s = s.optimal, type = "response")
    if (family == "poisson") 
        pearson.resid = sum(w * (yy - fitted)^2/fitted)
    if (family == "binomial") 
        pearson.resid = sum(w * (yy - fitted)^2/(fitted * (1 - 
            fitted)))
    (abs(pearson.resid - target))^2
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ glm }
\keyword{ geographically weighted regression }
\keyword{ gwr }
\keyword{ variable selection }

