library(mgcv)
library(spgwr)

setwd("~/git/gwr/code")
source("cubic.r")

data(trees)

rg = range(trees$Girth)
rh = range(trees$Height)

trees$Girth = (trees$Girth - rg[1]) / (rg[2] - rg[1])
trees$Height = (trees$Height - rh[1]) / (rh[2] - rh[1])

amT = am.setup(list(g=trees$Girth, h=trees$Height))

tree.model = fit.am(trees$Volume, amT[['X']], amT[['S']], c(0.01, 5000))