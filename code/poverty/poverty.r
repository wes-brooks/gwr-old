library(gwselect)
registerCores(n=7)

library(gridExtra)

#Import poverty data
pov = read.csv("~/git/gwr/data/poverty/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
pov$X09pop = as.numeric(gsub(",", "", as.character(pov$X09pop)))
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty',
    logitfampov='logit( proportion families in poverty)', pfire='proportion financial, insurance, real estate')

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

#Convert pov2 from a list to a data frame:
pov2 = data.frame(pov2)

#Correct the Y2K bug
pov2 = within(pov2, year <- as.numeric(as.character(year)) + 1900)
pov2 = within(pov2, year <- ifelse(year<1960, year+100, year))

pops = read.table("~/git/gwr/data/poverty/historicalpops.txt", headers=TRUE)

plots = list()
model = list()
bw = list()

plots.logistic = list()
model.logistic = list()
bw.logistic = list()

for (year in c(1960, 1970, 1980, 1990, 2000, 2006)) {
    #Use the lasso for GWR models of poverty with 2006 data:
    df = pov2[pov2$year==year,]
    df = merge(x=df, y=pops, by.x="fips", by.y="FIPS")
    weights = df[,paste("pop", year, sep="")]    

    #Define which variables we'll use as predictors of poverty:
    predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
    f = as.formula(paste("logitindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))
    bw[[as.character(year)]] = gwlars.sel(formula=f, data=df, coords=df[,c('x','y')], weights=weights, longlat=TRUE, gweight=bisquare, mode.select='AIC', method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)
    model[[as.character(year)]] = gwlars(formula=f, data=df, N=1, coords=df[,c('x','y')], weights=weights, longlat=TRUE, gweight=bisquare, bw=bw[[as.character(year)]], mode.select='AIC', s=NULL, method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)

    f = as.formula(paste("pindpov ~ -1 + ", paste(predictors, collapse="+"), sep=""))
    bw.logistic[[as.character(year)]] = gwglmnet.sel(formula=f, data=df, family='binomial', coords=df[,c('x','y')], weights=weights, longlat=TRUE, gweight=bisquare, mode.select='AIC', method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)
    model.logistic[[as.character(year)]] = gwglmnet(formula=f, data=df, N=1, family='binomial', coords=df[,c('x','y')], weights=weights, longlat=TRUE, gweight=bisquare, bw=bw.logistic[[as.character(year)]], mode.select='AIC', s=NULL, method="knn", tol=0.001, parallel=TRUE, precondition=FALSE, adapt=TRUE, verbose=FALSE)

    #Make a map using polygons provided to the plot function.              
    #Put the county names into a form that can be matched.
    county = map_data('county')
    
    plots[[as.character(year)]] = list()
    for (v in predictors) {
        plots[[as.character(year)]][[v]] = plot.gwselect(model[[as.character(year)]], var=v, polygons=county, title=v) + opts(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
    }

    plots.logistic[[as.character(year)]] = list()
    for (v in predictors) {
        plots.logistic[[as.character(year)]][[v]] = plot.gwselect(model.logistic[[as.character(year)]], var=v, polygons=county, title=v) + opts(plot.margin=unit(c(0,0,0,1), "cm")) + scale_x_continuous('') + scale_y_continuous('')
    }
    
    pp = plots[[as.character(year)]]
    pdf(paste('figures/poverty/', year, '.linear.coefficients.pdf', sep=''), width=8, height=16)
    grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], pp[['pwh']], pp[['pblk']], pp[['phisp']], pp[['metro']], ncol=2)
    dev.off()

    pp = plots.logistic[[as.character(year)]]
    pdf(paste('figures/poverty/', year, '.logistic.coefficients.pdf', sep=''), width=8, height=16)
    grid.arrange(pp[['pag']], pp[['pex']], pp[['pman']], pp[['potprof']], pp[['pfire']], pp[['pserve']], pp[['pwh']], pp[['pblk']], pp[['phisp']], pp[['metro']], ncol=2)
    dev.off()
}