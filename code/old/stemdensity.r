setwd("/afs/cs.wisc.edu/u/b/r/brook/private/research/paleon/stemdensity/stemdensity/code")

#load necessary R libraries
library(mgcv)
library(reshape)
library(plotrix)

#Set the working directory and import the code we'll use to make heatmaps.
source("matplot.r")
models = list()

#Import the stem density, count, and standard deviation data
stemdensity_raw <- read.csv("../data/stem.density.wi.csv", header=T)
stemdensity_stdev <- read.csv("../data/stem.density.sd.wi.csv", header=T)
counts <- read.csv("../data/count.wi.csv", header=T)

#Total the counts for individual taxa
stemdensity_raw$tot = apply(stemdensity_raw[,3:30], 1, sum)
counts$tot = apply(counts[,3:30], 1, sum)

#Remove spurious columns
stemdensity_raw <- within(stemdensity_raw, remove(X))
stemdensity_stdev <- within(stemdensity_stdev, remove(X))
counts <- within(counts, remove(X))

#Identify rows with no data
spp = names(stemdensity_raw)[!(names(stemdensity_raw) %in% c("Latitude", "Longitude"))]
empty = which(apply(stemdensity_raw[,spp], 1, function(x) sum(!is.na(x)))==0)

#Remove the rows with no data
stemdensity = stemdensity_raw[-empty,]
stemdensity_stdev = stemdensity_stdev[-empty,]
counts = counts[-empty,]

#Find the variance of the aggregate stem density by summing the variance of individual taxa
varstemdens = apply(stemdensity_stdev, 1, function(x) {sum(ifelse(x>0, x**2, 0))})

#Remove rows for which we have no variance data
indx = as.numeric(which(!is.na(varstemdens)))
varstemdens = varstemdens[indx]
stemdensity = stemdensity[indx,]

#Normalize the weights to average 1: this greatly speeds up the GAM-fitting
l = length(varstemdens)
tot = sum(1 / varstemdens)
w = l / tot
stemdensity$weight = w / varstemdens

#Make GAM model of the weighted stem density.
models[["gamma"]] = gam(tot ~ s(Latitude, Longitude, k=250), data=stemdensity, family=Gamma(link=log), weights=weight, lambda=1.4)

#Set the directory to store plots
plot_dir = "../figures/"
model = "gamma"
genus = "tot"
xx = range(stemdensity[,genus], na.rm=TRUE)


#Put the observed stem density and the weight into matrices that we can plot as heatmaps
#Create the matrix for the observed stem density (filled by default with NAs):
loc = with(stemdensity, list(lat=unique(Latitude), long=unique(Longitude)))
trees = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(trees) <- sort(unique(stemdensity$Latitude), decreasing=F)
colnames(trees) <- sort(unique(stemdensity$Longitude), decreasing=F)

#Create the matrix for the weights (filled by default with NAs):
loc = with(stemdensity, list(lat=unique(Latitude), long=unique(Longitude)))
weight_mat = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(weight_mat) <- sort(unique(stemdensity$Latitude), decreasing=F)
colnames(weight_mat) <- sort(unique(stemdensity$Longitude), decreasing=F)

#Put the observed stem densities and weights into their lat-long matrices
for(row in 1:dim(stemdensity)[1])
{
    trees[as.character(stemdensity[row, "Latitude"]), as.character(stemdensity[row, "Longitude"])] = ifelse((!is.nan(stemdensity[row, genus]) & stemdensity[row, genus]>0), stemdensity[row, genus], NA)

    weight_mat[as.character(stemdensity[row, "Latitude"]), as.character(stemdensity[row, "Longitude"])] = ifelse((!is.nan(stemdensity[row, genus]) & stemdensity[row, "weight"]>0), stemdensity[row, "weight"], NA)
}

#Put the model's fitted values into a matrix that we can plot as a heatmap (indexed by Lat, Long):
#First, create a data frame with Lat, Long, and fitted stem density
rows = which(!is.na(stemdensity[, genus]))
fits = cbind(stemdensity[rows, c("Latitude", "Longitude")], fitted=models[[model]]$fit)

#Now create a matrix of the appropriate dimensions (filled with NAs by default)
fitted = matrix(NA, nrow=length(loc[['lat']]), ncol=length(loc[['long']]))
rownames(fitted) <- sort(unique(fits$Latitude), decreasing=F)
colnames(fitted) <- sort(unique(fits$Longitude), decreasing=F)

#Put the fittes stem densities into the lat-long matrix
for(row in 1:dim(stemdensity)[1])
    fitted[as.character(fits$Latitude[row]), as.character(fits$Longitude[row])] = fits$fitted[row]                  


#Make several plots of the model and the data:
pdf(paste(plot_dir, "model_log_heatmap_", model, ".pdf", sep=""))
par(bty='n')
matplot(log(fitted), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=log(xx))
title(paste("log(expected[stem density]) from ", model, " model", sep=""))
dev.off()

pdf(paste(plot_dir, "model_heatmap_", model, ".pdf", sep=""))
par(bty='n')
matplot(fitted, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=xx)
title(paste("expected[stem density] from ", model, " model", sep=""))
dev.off()

pdf(paste(plot_dir, "model_stdev_heatmap_", model, ".pdf", sep=""))
par(bty='n')
matplot(sqrt(fitted**2 * models[[model]]$scale), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
title(paste("stdev(stemdensity) from ", model, " model", sep=""))
dev.off()

pdf(paste(plot_dir, "model_log_stdev_heatmap_", model, ".pdf", sep=""))
par(bty='n')
matplot(log(sqrt(fitted**2 * models[[model]]$scale)), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
title(paste("log(stdev[stem density]) from ", model, " model", sep=""))
dev.off()

pdf(paste(plot_dir, "model_Elog_heatmap_", model, ".pdf", sep=""))
par(bty='n')
matplot(log(fitted*models[[model]]$scale) + digamma(1/models[[model]]$scale), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=log(xx))
title(paste("expected[log(stem density)] from ", model, " model\n(variance is constant on this scale)", sep=""))
dev.off()

pdf(paste(plot_dir, "diagnostics", ".pdf", sep=""))
gam.check(models[[model]])
dev.off()                    

#plot the heatmap of stemdensity
pdf(paste(plot_dir, "data_heatmap.pdf", sep=""))
par(bty='n')
matplot(trees, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=xx)
title("measured stem density")
dev.off()

#plot the heatmap of stemdensity
pdf(paste(plot_dir, "log_data_heatmap.pdf", sep=""))
par(bty='n')
matplot(log(trees), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F, xrange=log(xx))
title("log of measured stem density")
dev.off()

#plot the heatmap of the weights
pdf(paste(plot_dir, "log_weight_heatmap.pdf", sep=""))
par(bty='n')
matplot(log(weight_mat), c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
title("log of weight used in fitting the GAM models")
dev.off()

#plot the heatmap of the weights
pdf(paste(plot_dir, "weight_heatmap.pdf", sep=""))
par(bty='n')
matplot(weight_mat, c(1,1), c(1,0), c(1,0), border=NA, show.legend=T, yrev=F, axes=F, ann=F)
title("weight used in fitting the GAM models")
dev.off()

pdf(paste(plot_dir, "density_histogram.pdf", sep=""))
hist(stemdensity[,genus], breaks=30)
dev.off()
