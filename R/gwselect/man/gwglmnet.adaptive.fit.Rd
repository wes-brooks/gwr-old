\name{gwglmnet.adaptive.fit}
\alias{gwglmnet.adaptive.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use the adaptive LASSO to fit a GLM in the GWR setting.
}
\description{
Use the adaptive LASSO to fit a GLM in the GWR setting.
}
\usage{
gwglmnet.adaptive.fit(x, y, coords, weight.matrix, s, verbose, family, prior.weights)
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
  \item{weight.matrix}{
%%     ~~Describe \code{weight.matrix} here~~
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
function (x, y, coords, weight.matrix, s, verbose, family, prior.weights) 
{
    gwglmnet.object = list()
    coords.unique = unique(coords)
    model = list()
    s.optimal = vector()
    adapt.normx = list()
    adapt.scale = list()
    cv.error = list()
    coef.scale = list()
    glm.step = list()
    for (i in 1:dim(coords.unique)[1]) {
        colocated = which(coords[, 1] == coords.unique[i, 1] & 
            coords[, 2] == coords.unique[i, 2])
        loow = weight.matrix[i, -colocated]
        if (sum(loow) == 0) {
            return(list(cv.error = Inf))
        }
        prior.loow = prior.weights[-colocated]
        reps = length(colocated)
        w <- prior.loow * loow
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
        else {
            m <- ncol(xx)
            n <- nrow(xx)
            one <- rep(1, n)
            meanx <- drop(one \%*\% xx)/n
            x.centered <- scale(xx, meanx, FALSE)
            normx <- sqrt(drop(one \%*\% (x.centered^2)))
            adapt.normx[[i]] = normx
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
            out.glm = try(glm(yy ~ xs, family = family, weights = w))
            if (class(out.glm) == "try-error") {
                cat(paste("Had to use the last glm for location ", 
                  i, "\n", sep = ""))
                glm.step[[i]] = out.glm = glm.step[[i - 1]]
            }
            else {
                glm.step[[i]] = out.glm
            }
            beta.glm = out.glm$coeff[2:(m + 1)]
            adapt.weight = abs(beta.glm)
            adapt.scale[[i]] = adapt.weight
            for (k in 1:dim(x.centered)[2]) {
                if (!is.na(adapt.weight[k])) {
                  xs[, k] = xs[, k] * adapt.weight[k]
                }
                else {
                  xs[, k] = rep(0, dim(xs)[1])
                  adapt.weight[k] = 0
                }
            }
            coef.scale[[i]] = adapt.weight/normx
            names(coef.scale[[i]]) = sapply(strsplit(names(coef.scale[[i]]), 
                "xs"), function(x) {
                x[2]
            })
            if (sum(coef.scale[[i]]) < 1e-10) {
                if (verbose) {
                  cat(paste("opted for the intercept-only model at location: ", 
                    i, "\n", sep = ""))
                }
                model[[i]] = NULL
                predictions = rep(coef(out.glm)[["(Intercept)"]], 
                  length(colocated))
                cv.error[[i]] = abs(matrix(predictions - matrix(y[colocated], 
                  nrow = reps, ncol = length(s))))
                s.optimal = c(s.optimal, max(s))
            }
            else {
                if (family == "binomial") {
                  model[[i]] = glmnet(x = xs, y = cbind(1 - yy, 
                    yy), lambda = s, family = family, weights = w)
                }
                else {
                  model[[i]] = glmnet(x = xs, y = yy, lambda = s, 
                    family = family, weights = w)
                }
                predictions = predict(model[[i]], newx = scale(matrix(x[colocated, 
                  ], nrow = reps, ncol = dim(xx)[2]), center = meanx, 
                  scale = normx/adapt.weight), type = "response", 
                  s = s)
                cv.error[[i]] = colSums(abs(matrix(predictions - 
                  matrix(y[colocated], nrow = reps, ncol = length(s)), 
                  nrow = reps, ncol = length(s))))
                s.optimal = c(s.optimal, s[which.min(cv.error[[i]])])
            }
        }
        if (verbose) {
            cat(paste(i, "\n", sep = ""))
        }
    }
    gwglmnet.object[["coef.scale"]] = coef.scale
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
