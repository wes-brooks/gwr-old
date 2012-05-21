library(spgwr)
library(lars)
library(maps)


#extract reference data
mapcounties <- map_data("county")
mapstates <- map_data("state")


#Import the plotting functions:
setwd("~/git/gwr/code")
source("utils.r")

#Import poverty data
pov = read.csv("~/git/gwr/data/upMidWestpov_Iowa_cluster_names.csv", header=TRUE)
years = c('60', '70', '80', '90', '00', '06')
column.map = list(pindpov='proportion individuals in poverty', 
    logitindpov='logit( proportion individuals in poverty )', pag='pag', pex='pex', pman='pman', 
    pserve='pserve', potprof='potprof', pwh='proportion white', pblk='proportion black', pind='pind',
    phisp='proportion hispanic', metro='metro', pfampov='proportion families in poverty',
    logitfampov='logit( proportion families in poverty)', pfire='pfire')

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




#Limit the data to just 100 points for now for computational reasons
pov2006 = pov2[pov2$year==2006,]
df = pov2006






#Use this trick to compute the matrix of distances very quickly
bw=3
n = dim(pov2006)[1]
D1 = matrix(rep(pov2006$x,n), n,n)
D2 = matrix(rep(pov2006$y,n), n,n)
D = sqrt((D1-t(D1))**2 + (D2-t(D2))**2)
w = bisquare(D, bw=bw)







#Define which variables we'll use as predictors of poverty:
predictors = c('pag', 'pex', 'pman', 'pserve', 'pfire', 'potprof', 'pwh', 'pblk', 'phisp', 'metro')
f = as.formula(paste("logitindpov ~ ", paste(predictors, collapse="+"), sep=""))

#Make a new variable with the name of each predictor:
for (col in predictors) {
    assign(col, vector())
}


model.data = df[,predictors]
model.data[['logitindpov']] = df[['logitindpov']]





cv_error = data.frame()
w.lasso.geo = list()
coefs = list()
ss = seq(0, 1, length.out=100)
lambda = seq(0, 2, length.out=2000)
l = vector()



for(i in 1:dim(df)[1]) {
    loow = w[i,-i]

    model = lm(f, data=model.data[-i,], weights=loow)
    
    w.eig <- eigen(diag(loow))
    w.sqrt <- w.eig$vectors %*% diag(sqrt(w.eig$values)) %*% solve(w.eig$vectors)
    w.lasso.geo[[i]] = lars(x=w.sqrt %*% as.matrix(df[-i,predictors]), y=as.matrix(df$logitindpov[-i]))
    
    for (col in predictors) {
        coefs[[col]] = c(coefs[[col]], model$coef[[col]])
    }
    
    l = c(l, which.min(abs(predict(w.lasso.geo[[i]], newx=model.data[i,-11], s=lambda, type='fit', mode='lambda')[['fit']] - model.data[i,11]))/1000)
    print(i)
}



#Prepare something for plotting:
output = vector()
var = 'pserve'
df.temp = df
for (i in 1:dim(df.temp)[1]) {
    output = c(output, coef.lars(w.lasso.geo[[i]], newx=model.data[i,-11], mode='lambda', s=1000*l[i])[[var]])
}




df.temp$output = output
df.temp$lambda = l




#define color buckets
colors = c("#980043", "#C994C7", "#D4B9DA", "#DD1C77", "#DF65B0", "#F1EEF6")

bin = vector()
breaks = c(0.1, 0.20, 0.3, 0.4, 0.5, 2.00)
for (i in 1:dim(df.temp)[1]) {
    bin = c(bin, min(which(breaks >= df.temp$lambda[i])))
}
df.temp$color = colors[bin]



df.temp$county = tolower(as.character(df.temp$COUNTY))
for (i in 1:dim(df.temp)[1]) {
    county = gsub("['-. ]", '', df.temp$county[i])
    df.temp$county[i] = paste(county, tolower(df$STATE[i]), sep=',')
}


#extract reference data
mapcounties <- map_data("county")
mapstates <- map_data("state")

#limit our view to the midwest:
midweststates = mapstates[tolower(mapstates$region) %in% tolower(df.temp$STATE),]
midwestcounties = mapcounties[tolower(mapcounties$region) %in% tolower(df.temp$STATE),]

#merge data with ggplot county coordinates
midwestcounties$county <- with(midwestcounties , paste(gsub("['-. ]", '', subregion), region, sep = ","))
mergedata <- merge(midwestcounties, df.temp, by.x = "county", by.y = "county")
mergedata <- mergedata[order(mergedata$group, mergedata$order),]

#draw map
map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=output))
map <- map + scale_fill_gradient(low='white', high='red', limits=range(output)) + #scale_fill_brewer(palette="PuRd") +
    coord_map(project="globular") +
    opts(legend.position = "none")

map <- map + opts(panel.background = theme_rect(fill='green', colour='red'))

#add state borders
map <- map + geom_path(data = midweststates, colour = "white", size = .75)

#add county borders
map <- map + geom_path(data = midwestcounties, colour = "white", size = .5, alpha = .1)
map

