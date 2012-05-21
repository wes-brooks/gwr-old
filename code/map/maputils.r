unemp = read.csv("../data/unemployment09.csv", header=FALSE)
abbrevs = list('alabama'='al', 'alaska'='ak', 'arizona'='az', 'arkansas'='ar', 'california'='ca', 'colorado'='co', 'connecticut'='ct', 'delaware'='de', 'florida'='fl', 'georgia'='ga', 'hawaii'='hi', 'idaho'='id', 'illinois'='il', 'indiana'='in', 'iowa'='ia', 'kansas'='ks', 'kentucky'='ky', 'louisiana'='la', 'maine'='me', 'maryland'='md', 'massachusetts'='ma', 'michigan'='mi', 'minnesota'='mn', 'mississippi'='ms', 'missouri'='mo', 'montana'='mt', 'nebraska'='ne', 'nevada'='nv', 'new hampshire'='nh', 'new jersey'='nj', 'new mexico'='nm', 'new york'='ny', 'north carolina'='nc', 'north dakota'='nd', 'ohio'='oh', 'oklahoma'='ok', 'oregon'='or', 'pennsylvania'='pa', 'rhode island'='ri', 'south carolina'='sc', 'south dakota'='sd', 'tennessee'='tn', 'texas'='tx', 'utah'='ut', 'vermont'='vt', 'virginia'='va', 'washington'='wa', 'west virginia'='wv', 'wisconsin'='wi', 'wyoming'='wy', 'puerto rico'='pr', 'columbia'='dc')
states = as.vector(names(abbrevs), mode='list')
names(states) = as.vector(unlist(abbrevs))
descriptors = c( 'county/city', 'county/town', 'borough/municipality', 'borough/city', 'county', 'parish', 'precinct', 'ward', 'borough', 'city', 'municipio', 'census area')

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
locality =  strsplit(tolower(unemp$county), ', ')
for (i in 1:length(locality)) {
    for (j in 1:length(descriptors)) {
        locality[[i]][1] = sub(descriptors[j], '', locality[[i]][1])
    }
    locality[[i]][1] = gsub("(miami)?['-. ]", '', locality[[i]][1])
    unemp$county[i] = paste(locality[[i]][1], states[[locality[[i]][2]]], sep=',')
}
unemp$stateid = as.factor(unemp$stateid)
unemp$countyid = as.factor(unemp$countyid)

#unemp$colorBins = colors[ceiling(length(colors)*rank(-unemp$unemployment)/dim(unemp)[1])]


#define color buckets
colors = c("#980043", "#C994C7", "#D4B9DA", "#DD1C77", "#DF65B0", "#F1EEF6")

bin = vector()
breaks = c(2, 4, 6, 8, 10, 100)
for (i in 1:dim(unemp)[1]) {
    bin = c(bin, min(which(breaks >= unemp$unemployment[i])))
}
unemp$color = colors[bin]



#extract reference data
mapcounties <- map_data("county")
mapstates <- map_data("state")

#merge data with ggplot county coordinates
mapcounties$county <- with(mapcounties , paste(gsub("(miami)?['-. ]", '', subregion), region, sep = ","))
mergedata <- merge(mapcounties, unemp, by.x = "county", by.y = "county")
mergedata <- mergedata[order(mergedata$group, mergedata$order),]

#draw map
map <- ggplot(mergedata, aes(long,lat,group=group)) + geom_polygon(aes(fill=color))
map <- map + scale_fill_brewer(palette="PuRd") +
    coord_map(project="globular") +
    opts(legend.position = "none")

map <- map + opts(panel.background = theme_rect(fill='green', colour='red'))

#add state borders
map <- map + geom_path(data = mapstates, colour = "white", size = .75)

#add county borders
map <- map + geom_path(data = mapcounties, colour = "white", size = .5, alpha = .1)
map