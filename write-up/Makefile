SHELL = /bin/sh

all: template.pdf

template.pdf: *.tex *.bib */*.tex
	pdflatex template.tex
	bibtex template 
	touch template.tex
	pdflatex template.tex
	pdflatex template.tex

short: *.tex *.bib */*.tex
	pdflatex template.tex

clean:
	if (rm -f *.log *.aux *.bbl *.blg *.toc *.dvi *.ps *.brf *.lbl *pdf) then :; fi

