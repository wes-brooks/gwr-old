.PHONY : tables
.PHONY : extract

#Default variable values (in case none were specified on the command line):
outdir ?= ~/git/gwr/output/simulation/
root ?= ~/misc/gwr/
cluster ?= 96


#Import the saved Rdata and turn it into simulation tables.
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
	  --header=selection-header.tex --sizex=16cm
	python3 format-tables.py --infile $(outdir)/../prelim-talk/msex-step.tex --outfile $(outdir)/../prelim-talk/msex-step.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/msey-step.tex --outfile $(outdir)/../prelim-talk/msey-step.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/bx-step.tex --outfile $(outdir)/../prelim-talk/bx-step.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/by-step.tex --outfile $(outdir)/../prelim-talk/by-step.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/varx-step.tex --outfile $(outdir)/../prelim-talk/varx-step.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/vary-step.tex --outfile $(outdir)/../prelim-talk/vary-step.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/msex-gradient.tex --outfile $(outdir)/../prelim-talk/msex-gradient.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/msey-gradient.tex --outfile $(outdir)/../prelim-talk/msey-gradient.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/bx-gradient.tex --outfile $(outdir)/../prelim-talk/bx-gradient.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/by-gradient.tex --outfile $(outdir)/../prelim-talk/by-gradient.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/varx-gradient.tex --outfile $(outdir)/../prelim-talk/varx-gradient.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/vary-gradient.tex --outfile $(outdir)/../prelim-talk/vary-gradient.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/msex-parabola.tex --outfile $(outdir)/../prelim-talk/msex-parabola.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/msey-parabola.tex --outfile $(outdir)/../prelim-talk/msey-parabola.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/bx-parabola.tex --outfile $(outdir)/../prelim-talk/bx-parabola.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/by-parabola.tex --outfile $(outdir)/../prelim-talk/by-parabola.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/varx-parabola.tex --outfile $(outdir)/../prelim-talk/varx-parabola.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/vary-parabola.tex --outfile $(outdir)/../prelim-talk/vary-parabola.tex --scale=0.55
	python3 format-tables.py --infile $(outdir)/../prelim-talk/selection.tex \
	  --outfile $(outdir)/../prelim-talk/selection.tex --tabcolsep=4pt \
	  --header=selection-header.tex --scale=0.6


#extract the tarred-and-gzipped simulation results, import them to R, and save the resulting Rdata.
extract: 
	mkdir -p $(root)/output
	sh extract.sh $(root) $(cluster)
	python3 cleanup-results.py $(root) $(cluster)
	Rscript summarize.r $(root) $(cluster)

