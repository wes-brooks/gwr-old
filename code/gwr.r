library(spgwr)
data(columbus)

pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
heads = c('pindpov', 'logitindpov', 'pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'pind', 'phisp', 'metro', 'pfampov', 'logitfampov')
years = c('60', '70', '80', '90', '00', '06')

pov2 = list()

for (column.name in heads) {
    col = vector()

    for (year in years) {
        if (paste(column.name, year, sep="") %in% names(pov)) {
            indx = which(names(pov)==paste(column.name, year, sep=""))
            col = c(col, pov[,indx])
        }
        else {
            col = c(col, rep(NA, dim(pov)[1]))
        }
    }
    pov2[[column.name]] = col
}

#Find the columns we haven't yet matched:
"%w/o%" <- function(x, y) x[!x %in% y]
missed = names(pov) %w/o% outer(heads, years, FUN=function(x, y) {paste(x, y, sep="")})

for (column.name in missed) {
    col = rep(pov[,column.name], length(years))
    pov2[[column.name]] = col
}

pov2[['year']] = vector()
for (year in years) {
    pov2[['year']] = c(pov2[['year']], rep(year, dim(pov)[1]))
}

pov2 = data.frame(pov2)

pov2 = within(pov2, year <- as.numeric(as.character(year)) + 1900)
pov2 = within(pov2, year <- ifelse(year<1960, year+100, year))

bw = gwr.sel(income~1, data=columbus, coords=cbind(x,y), adapt=FALSE, gweight=gwr.bisquare)
gwr.model1 = gwr(income~1, data=columbus, coords=cbind(x,y), bandwidth=bw, gweight=gwr.bisquare, hatmatrix=TRUE)

