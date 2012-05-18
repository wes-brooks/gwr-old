library(ggplot2)


#load data
stateAbbr rownames(stateAbbr)unemp_data

#get county names in correct format
countyNames counties statesstates

#concatenate states and counties
unemp_data$counties

#parse out county titles & specifics
unemp_data$counties unemp_data$counties unemp_data$countiesunemp_data$counties

#define color buckets
colors = c("#F1EEF6", "#D4B9DA", "#C994C7", "#DF65B0", "#DD1C77", "#980043")
unemp_data$colorBuckets



#extract reference data
mapcounties <- map_data("county")
mapstates <- map_data("state")

#merge data with ggplot county coordinates
mapcounties$county <- with(mapcounties , paste(subregion, region, sep = ","))
mergedata <- merge(mapcounties, unemp, by.x = "county", by.y = "county")
mergedata <- mergedata[order(mergedata$order),]

#draw map
map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=colorBuckets))
map <- map + scale_fill_brewer(palette="PuRd") +
    coord_map(project="globular") +
    opts(legend.position = "none")

#add state borders
map <- map + geom_path(data = mapstates, colour = "white", size = .75)

#add county borders
map <- map + geom_path(data = mapcounties, colour = "white", size = .5, alpha = .1)
map


unemp = read.csv("../data/unemployment09.csv", header=FALSE)

pop1 = vector()
pop2 = vector()
pop3 = vector()
for (i in 1:dim(unemp)[1]) {
    pop1 = c(pop1, as.numeric(paste(strsplit(sub("[[:blank:]]+$", "", as.character(unemp$V6[i])), ',')[[1]], collapse='')))
    pop2 = c(pop2, as.numeric(paste(strsplit(sub("[[:blank:]]+$", "", as.character(unemp$V7[i])), ',')[[1]], collapse='')))
    pop3 = c(pop3, as.numeric(paste(strsplit(sub("[[:blank:]]+$", "", as.character(unemp$V8[i])), ',')[[1]], collapse='')))
}
unemp$V6 = pop1
unemp$V7 = pop2
unemp$V8 = pop3

names(unemp) = c("LAUS", 'stateid', 'countyid', 'county', 'year', '2009 population', '2000 population', 'population change', 'unemployment')
unemp$year = as.numeric(unemp$year) 
unemp$county = as.character(unemp$county)
unemp$county = sapply(1:length(unemp$county), FUN=function(x) {tolower(paste(strsplit(unemp$county[x], " ")[[1]][1], tail(strsplit(unemp$county[x], " ")[[1]],1), sep=', '))})
unemp$stateid = as.factor(unemp$stateid)
unemp$countyid = as.factor(unemp$countyid)



unemp$colorBins = colors[ceiling(length(colors)*order(unemp$unemployment)/dim(unemp)[1])]

map <- ggplot() + layer(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=colorBins))
map <- map + scale_fill_brewer(palette="PuRd") +
    coord_map(project="globular") +
    opts(legend.position = "none")