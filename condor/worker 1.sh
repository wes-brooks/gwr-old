mkdir output
#mkdir logs
Rscript code/simulations.r $*

tar -zcvf output/output$1-$2.tgz2 output/*.csv
rm output/*.csv
