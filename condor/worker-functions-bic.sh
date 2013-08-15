if [ -f /usr/bin/R ]
then
	mkdir output

	if $(uname -m |grep '64')
	then
		if $(R --version | grep 'R version 2')
		then
			Rscript code/simulations-functions-bic-x86_64.r $*
		else
			Rscript code/simulations-functions-bic-3-x86_64.r $*
		fi
	else
		if $(R --version | grep 'R version 2')
		then
			Rscript code/simulations-functions-bic-i386.r $*
		else
			Rscript code/simulations-functions-bic-3-i386.r $*
		fi
	fi

	tar -zcvf output/output-$1-$2.tgz2 output/*.csv
	rm output/*.csv

	exit 0
else 
         exit 1
fi
