library(gwselect, lib.loc=c('R', 'R/x86_64-redhat-linux-gnu-library/2.15'))
args = commandArgs(trailingOnly=TRUE)
write(paste('hello: ', paste(args, collapse=','), '\n', sep=''), file='output/hello.txt', append=TRUE)
