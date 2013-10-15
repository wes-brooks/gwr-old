poverty-plots: plots.r poverty.r
	Rscript -e "require(RCurl); \
	  require(devtools); \
	  install_github('wesesque/brooks'); \
	  brooks::source_https('https://raw.github.com/wesesque/gwr/master/code/poverty/poverty.r'); \
	  brooks::source_https('https://raw.github.com/wesesque/gwr/master/code/poverty/plots.r');"
