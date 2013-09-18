.PHONY : tables
.PHONY : extract

tables:
	Rscript summarize.r $(root) $(cluster)
	Rscript summary-tables.r $(root) $(cluster) $(outdir)
	python3 format-tables.py --infile $(outdir)/msex.tex --outfile $(outdir)/msex.tex
	python3 format-tables.py --infile $(outdir)/msey.tex --outfile $(outdir)/msey.tex
	python3 format-tables.py --infile $(outdir)/bx.tex --outfile $(outdir)/bx.tex
	python3 format-tables.py --infile $(outdir)/by.tex --outfile $(outdir)/by.tex
	python3 format-tables.py --infile $(outdir)/varx.tex --outfile $(outdir)/varx.tex
	python3 format-tables.py --infile $(outdir)/vary.tex --outfile $(outdir)/vary.tex
	python3 format-tables.py --infile $(outdir)/selection.tex --outfile $(outdir)/selection.tex --boxsize=22cm


extract: 
	mkdir -p $(root)/output
	sh extract.sh $(root) $(cluster)
	python3 cleanup-results.py $(root) $(cluster)
