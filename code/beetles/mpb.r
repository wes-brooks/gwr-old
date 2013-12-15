require(RCurl)
require(devtools)
install_github("gwselect", "wesesque")
require(gwselect)
#registerCores(n=3)

load_https <- function(url, sep=',', header=TRUE, row.names=NULL, ...) {
  # Import the data:
  read.table(text = getURL(url,
    followlocation=TRUE, cainfo=system.file("CurlSSL", "cacert.pem", package="RCurl")),
    sep=sep, header=header, row.names=row.names, ...)
}

mpb = load_https("https://raw.github.com/wesesque/gwr/master/data/SouthernPineBeetle/Code-Andy/mpb.csv")

#Remove rows with NAs:
n = nrow(mpb)
indx = which(is.na(mpb))
na.rows = (indx-1) %% n + 1
mpb = mpb[-na.rows,]

predictors = c('meanelevation', 'warm', 'Tmin', 'Tmean', 'Tmax', 'cold', 'precip', 'dd', 'ddegg')
f = as.formula(paste("nifestations ~ -1 + ", paste(predictors, collapse="+"), sep=""))
bw = gwglmnet.sel(formula=f, data=mpb, family='poisson', alpha=1, coords=mpb[,c('X','Y')], longlat=FALSE, mode.select="BIC", gweight=spherical, tol=1, s=NULL, method='dist', adapt=TRUE, parallel=FALSE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

sink("~/mpb.out.txt")
print(bw)
sink()

#model = gwglmnet.nen(nifestations~meanelevation+warm+Tmin+Tmean+Tmax+cold+precip+dd+ddegg, data=mpb, coords=mpb[,c('X','Y')], gweight=bisquare, s=seq(0,5,0.001), tol=10, bw=200000, type='pearson', family='poisson', parallel=TRUE, weights=weights)
