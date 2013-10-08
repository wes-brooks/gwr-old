openaccess.pdf: openaccess.tex
	xelatex openaccess

notes: openaccess_withnotes.pdf

all: openaccess.pdf notes

openaccess_withnotes.pdf: openaccess_withnotes.tex
	xelatex openaccess_withnotes
	pdfnup openaccess_withnotes.pdf --nup 1x2 --no-landscape --paper letterpaper --frame true --scale 0.9
	mv openaccess_withnotes-nup.pdf openaccess_withnotes.pdf

openaccess_withnotes.tex: openaccess.tex createVersionWithNotes.rb
	createVersionWithNotes.rb openaccess.tex openaccess_withnotes.tex

web: openaccess.pdf dropbox openaccess_withnotes.pdf
	scp openaccess*.pdf broman-2:public_html/presentations/

dropbox: ~/Dropbox/Talks/openaccess.pdf

~/Dropbox/Talks/openaccess.pdf: openaccess.pdf openaccess_withnotes.pdf
	cp openaccess*.pdf ~/Dropbox/Talks/
