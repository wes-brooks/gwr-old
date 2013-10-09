practice-talk.pdf: practice-talk.tex
	xelatex practice-talk

notes: practice-talk_withnotes.pdf

all: practice-talk.pdf notes

practice-talk_withnotes.pdf: practice-talk_withnotes.tex
	xelatex practice-talk_withnotes
	pdfnup practice-talk_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv practice-talk_withnotes-nup.pdf practice-talk_withnotes.pdf

practice-talk_withnotes.tex: practice-talk.tex createVersionWithNotes.rb
	./createVersionWithNotes.rb practice-talk.tex practice-talk_withnotes.tex

web: practice-talk.pdf dropbox practice-talk_withnotes.pdf
	scp practice-talk*.pdf wesley@web318.webfaction.com:/home/wesley/slideshows/

dropbox: ~/Dropbox/Talks/practice-talk.pdf

~/Dropbox/Talks/practice-talk.pdf: practice-talk.pdf practice-talk_withnotes.pdf
	cp practice-talk*.pdf ~/Dropbox/Talks/
