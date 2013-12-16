require(RCurl)
require(devtools)
install_github("gwselect", "wesesque")
require(gwselect)
registerCores(n=3)

load_https <- function(url, sep=',', header=TRUE, row.names=NULL, ...) {
  # Import the data:
  read.table(text = getURL(url,
    followlocation=TRUE, cainfo=system.file("CurlSSL", "cacert.pem", package="RCurl")),
    sep=sep, header=header, row.names=row.names, ...)
}

aphids = load_https("https://raw.github.com/wesesque/gwr/master/data/SoybeanAphids/data_8_covariates.csv")

#Remove rows with NAs:
n = nrow(aphids)
indx = which(is.na(aphids))
na.rows = (indx-1) %% n + 1
if (length(na.rows)>0) aphids = aphids[-na.rows,]

aphids$count = round(exp(aphids$log_count))

aphids06 = aphids[aphids$year==2006,]
n = nrow(aphids06)

predictors = c('Corn', 'soybeans', 'S_Grains', 'All_Grasslands', 'Forests_All', 'Water', 'Developed', 'Wetlands')
f = as.formula(paste("count ~ -1 + ", paste(predictors, collapse="+"), sep=""))
bw = gwglmnet.sel(formula=f, data=aphids06, family='poisson', alpha=1, coords=aphids06[,c('Longitude','Latitude'), longlat=TRUE, mode.select="BIC", gweight=spherical, tol=1, s=NULL, method='dist', adapt=TRUE, parallel=TRUE, interact=TRUE, verbose=TRUE, shrunk.fit=FALSE)

#model = gwglmnet.nen(nifestations~meanelevation+warm+Tmin+Tmean+Tmax+cold+precip+dd+ddegg, data=mpb, coords=mpb[,c('X','Y')], gweight=bisquare, s=seq(0,5,0.001), tol=10, bw=200000, type='pearson', family='poisson', parallel=TRUE, weights=weights)
