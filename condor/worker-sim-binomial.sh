if [ -f /unsup/R/bin/R ]
then
	mkdir output

    /unsup/R/bin/Rscript code/simulation-binomial.r $*

	tar -zcvf output/output-binomial-$1-$2.tgz2 output/*.csv
	rm output/*.csv

	exit 0
else 
     exit 1
fi
