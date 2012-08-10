landcover = read.csv("~/git/gwr/data/NorthernWisconsin/Data-Chongyang/landcover.csv", head=TRUE)
landcover$Cover1987 = as.numeric(ifelse(landcover$Cover1987=="Aspen-Paper Birch",1,0))
#source("~/git/gwr/code/ashland/convert_coords.r")
ohfive = read.csv("~/git/gwr/data/NorthernWisconsin/Data-Chongyang/1905.csv", head=TRUE)
landcover05 = cbind(landcover, ohfive)

source("~/git/gwr/code/utils.r")
source("~/git/gwr/code/gwglmnet.nen.sel.r")
#Remove rows with NAs:
n = nrow(landcover05)
indx = which(is.na(landcover05))
na.rows = (indx-1) %% n + 1
if (length(na.rows)>0) landcover05 = landcover05[-na.rows,]

n = nrow(landcover05)
weights = rep(1, n)

bw = gwglmnet.nen.sel(Cover1987~Own+Res+PolyNm+PolyPr+MxPolyPr+TotOwn+AvParcel, data=landcover05, coords=landcover05[,c('X','Y')], weights=weights, gweight=bisquare, s=NULL, tol=10, type='pearson', family='binomial', parallel=TRUE, adapt=TRUE)
