if [ -f /unsup/R/bin/R ]
then
	mkdir output

    /unsup/R/bin/Rscript code/simulation-gaussian.r $*

	tar -zcvf output/output-$1-$2.tgz2 output/*.csv
	rm output/*.csv

	exit 0
else 
     exit 1
fi
