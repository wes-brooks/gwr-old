.PHONY : tables
.PHONY : extract

tables:
	Rscript summarize.r $(root) $(cluster)
	Rscript summary-tables.r $(root) $(cluster) $(outdir)

extract: 
	mkdir -p $(root)/output
	sh extract.sh $(root) $(cluster)
	python3 cleanup-results.py $(root) $(cluster)
