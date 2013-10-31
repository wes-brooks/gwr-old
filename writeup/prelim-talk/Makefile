prelim-talk.pdf: prelim-talk.tex
	xelatex prelim-talk
	biber prelim-talk
	xelatex prelim-talk

notes: prelim-talk_withnotes.pdf

all: prelim-talk.pdf notes

prelim-talk_withnotes.pdf: prelim-talk_withnotes.tex
	xelatex prelim-talk_withnotes
	biber prelim-talk_withnotes
	xelatex prelim-talk_withnotes
	pdfnup prelim-talk_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv prelim-talk_withnotes-nup.pdf prelim-talk_withnotes.pdf

prelim-talk_withnotes.tex: prelim-talk.tex createVersionWithNotes.rb
	./createVersionWithNotes.rb prelim-talk.tex prelim-talk_withnotes.tex

web: prelim-talk.pdf dropbox prelim-talk_withnotes.pdf
	scp prelim-talk*.pdf wesley@web318.webfaction.com:/home/wesley/slideshows/

dropbox: ~/Dropbox/Talks/prelim-talk.pdf

~/Dropbox/Talks/prelim-talk.pdf: prelim-talk.pdf prelim-talk_withnotes.pdf
	cp prelim-talk*.pdf ~/Dropbox/Talks/
