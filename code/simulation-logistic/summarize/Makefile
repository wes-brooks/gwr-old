.PHONY : tables
.PHONY : extract

#Default variable values (in case none were specified on the command line):
outdir ?= ~/git/gwr/output/simulation-binomial
root ?= ~/gwr/simoutput-binomial


#Import the saved Rdata and turn it into simulation tables.
tables:
	Rscript summary-tables.r $(root) binomial $(outdir)
	python format-tables.py --infile $(outdir)/msex.tex --outfile $(outdir)/msex.tex --sizey=12cm
	python format-tables.py --infile $(outdir)/msey.tex --outfile $(outdir)/msey.tex --sizey=12cm
	python format-tables.py --infile $(outdir)/bx.tex --outfile $(outdir)/bx.tex --sizey=12cm
	python format-tables.py --infile $(outdir)/by.tex --outfile $(outdir)/by.tex --sizey=12cm
	python format-tables.py --infile $(outdir)/varx.tex --outfile $(outdir)/varx.tex --sizey=12cm
	python format-tables.py --infile $(outdir)/vary.tex --outfile $(outdir)/vary.tex --sizey=12cm
	python format-tables.py --infile $(outdir)/selection.tex \
	  --outfile $(outdir)/selection.tex --tabcolsep=4pt \
	  --header=selection-header.tex --sizex=16cm
	python format-tables.py --infile $(outdir)/msex-mini-step.tex --outfile $(outdir)/msex-mini-step.tex --scale=0.55
	python format-tables.py --infile $(outdir)/msey-mini-step.tex --outfile $(outdir)/msey-mini-step.tex --scale=0.55
	python format-tables.py --infile $(outdir)/bx-mini-step.tex --outfile $(outdir)/bx-mini-step.tex --scale=0.55
	python format-tables.py --infile $(outdir)/by-mini-step.tex --outfile $(outdir)/by-mini-step.tex --scale=0.55
	python format-tables.py --infile $(outdir)/varx-mini-step.tex --outfile $(outdir)/varx-mini-step.tex --scale=0.55
	python format-tables.py --infile $(outdir)/vary-mini-step.tex --outfile $(outdir)/vary-mini-step.tex --scale=0.55
	python format-tables.py --infile $(outdir)/msex-mini-gradient.tex --outfile $(outdir)/msex-mini-gradient.tex --scale=0.55
	python format-tables.py --infile $(outdir)/msey-mini-gradient.tex --outfile $(outdir)/msey-mini-gradient.tex --scale=0.55
	python format-tables.py --infile $(outdir)/bx-mini-gradient.tex --outfile $(outdir)/bx-mini-gradient.tex --scale=0.55
	python format-tables.py --infile $(outdir)/by-mini-gradient.tex --outfile $(outdir)/by-mini-gradient.tex --scale=0.55
	python format-tables.py --infile $(outdir)/varx-mini-gradient.tex --outfile $(outdir)/varx-mini-gradient.tex --scale=0.55
	python format-tables.py --infile $(outdir)/vary-mini-gradient.tex --outfile $(outdir)/vary-mini-gradient.tex --scale=0.55
	python format-tables.py --infile $(outdir)/msex-mini-parabola.tex --outfile $(outdir)/msex-mini-parabola.tex --scale=0.55
	python format-tables.py --infile $(outdir)/msey-mini-parabola.tex --outfile $(outdir)/msey-mini-parabola.tex --scale=0.55
	python format-tables.py --infile $(outdir)/bx-mini-parabola.tex --outfile $(outdir)/bx-mini-parabola.tex --scale=0.55
	python format-tables.py --infile $(outdir)/by-mini-parabola.tex --outfile $(outdir)/by-mini-parabola.tex --scale=0.55
	python format-tables.py --infile $(outdir)/varx-mini-parabola.tex --outfile $(outdir)/varx-mini-parabola.tex --scale=0.55
	python format-tables.py --infile $(outdir)/vary-mini-parabola.tex --outfile $(outdir)/vary-mini-parabola.tex --scale=0.55
	python format-tables.py --infile $(outdir)/selection.tex \
	  --outfile $(outdir)/selection.tex --tabcolsep=4pt \
	  --header=selection-header.tex --scale=0.6


#extract the tarred-and-gzipped simulation results, import them to R, and save the resulting Rdata.
extract: 
	mkdir -p $(root)/summary-binomial
	python cleanup-results.py $(root)
	Rscript summarize-binomial.r $(root) binomial

