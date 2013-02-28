if [ -f /usr/bin/R ]
then
	mkdir output
	
	Rscript code/simulations.r $*

	tar -zcvf output/output$1-$2.tgz2 output/*.csv
	rm output/*.csv
	exit 0
else 
         exit 1
fi