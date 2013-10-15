require(brooks)

#Import poverty data
pov = brooks::load_https("https://raw.github.com/wesesque/gwr/master/data/poverty/upMidWestpov_Iowa_cluster_names.csv", header=TRUE, sep=",", row.names=NULL)
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='Proportion of individuals in poverty', 
    logitindpov='logit of proportion of individuals in poverty',
    pag='Proportion working in agriculture',
    pex='Proportion working in mining',
    pman='Proportion working in manufacturing',
    pserve='Proportion working in services',
    potprof='Proportion working in other professions',
    pwh='Proportion white',
    pblk='Proportion black',
    pind='Proportion indiginous',
    phisp='Proportion hispanic',
    metro='metro',
    pfampov='Proportion of families in poverty',
    logitfampov='logit of proportion of families in poverty',
    pfire='Proportion working in\nfinance, insurance, or real estate'
)

#Process the poverty data so that each column appears only once and the year is added as a column.
pov2 = list()
for (column.name in names(column.map)) {
    col = vector()
    for (year in years) {
        if (paste(column.name, year, sep="") %in% names(pov)) {
            indx = which(names(pov)==paste(column.name, year, sep=""))
            col = c(col, pov[,indx])
        }
        else { col = c(col, rep(NA, dim(pov)[1])) }
    }
    pov2[[column.name]] = col
}

#Find the columns we haven't yet matched:
"%w/o%" <- function(x, y) x[!x %in% y]
missed = names(pov) %w/o% outer(names(column.map), years, FUN=function(x, y) {paste(x, y, sep="")})

for (column.name in missed) {
    col = rep(pov[,column.name], length(years))
    pov2[[column.name]] = col
}

#Add the year column to the pov2 data list.
pov2[['year']] = vector()
for (year in years) {
    pov2[['year']] = c(pov2[['year']], rep(year, dim(pov)[1]))
}

#Convert pov2 from a list to a data frame and correct the Y2K bug:
pov2 = data.frame(pov2)
pov2 = within(pov2, year <- as.numeric(as.character(year)) + 1900)
pov2 = within(pov2, year <- ifelse(year<1960, year+100, year))
