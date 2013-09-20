.PHONY : tables
.PHONY : extract

tables:
	Rscript summary-tables.r $(root) $(cluster) $(outdir)
	python3 format-tables.py --infile $(outdir)/msex.tex --outfile $(outdir)/msex.tex --sizey=12cm
	python3 format-tables.py --infile $(outdir)/msey.tex --outfile $(outdir)/msey.tex --sizey=12cm
	python3 format-tables.py --infile $(outdir)/bx.tex --outfile $(outdir)/bx.tex --sizey=12cm
	python3 format-tables.py --infile $(outdir)/by.tex --outfile $(outdir)/by.tex --sizey=12cm
	python3 format-tables.py --infile $(outdir)/varx.tex --outfile $(outdir)/varx.tex --sizey=12cm
	python3 format-tables.py --infile $(outdir)/vary.tex --outfile $(outdir)/vary.tex --sizey=12cm
	python3 format-tables.py --infile $(outdir)/selection.tex \
	  --outfile $(outdir)/selection.tex --tabcolsep=4pt \
	  --header=../../../output/simulation/selection-header.tex --sizex=16cm


extract: 
	mkdir -p $(root)/output
	sh extract.sh $(root) $(cluster)
	python3 cleanup-results.py $(root) $(cluster)
	Rscript summarize.r $(root) $(cluster)
