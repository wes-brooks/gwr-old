library(gwselect)

aphids = read.csv("~/git/gwr/data/SoybeanAphids/data_8 covariates.csv", head=TRUE)

#Remove rows with NAs:
n = nrow(aphids)
indx = which(is.na(aphids))
na.rows = (indx-1) %% n + 1
if (length(na.rows)>0) aphids = aphids[-na.rows,]

aphids$count = round(exp(aphids$log_count))

aphids06 = aphids[aphids$year==2006,]
n = nrow(aphids06)
weights = rep(1, n)

bw = gwglmnet.nen.sel(count~Corn + soybeans + S_Grains + All_Grasslands + Forests_All + Water + Developed + Wetlands, data=aphids06, coords=aphids06[,c('Longitude','Latitude')], weights=weights, gweight=bisquare, s=NULL, tol=10, type='pearson', family='poisson', parallel=TRUE, longlat=TRUE)

#model = gwglmnet.nen(nifestations~meanelevation+warm+Tmin+Tmean+Tmax+cold+precip+dd+ddegg, data=mpb, coords=mpb[,c('X','Y')], gweight=bisquare, s=seq(0,5,0.001), tol=10, bw=200000, type='pearson', family='poisson', parallel=TRUE, weights=weights)
