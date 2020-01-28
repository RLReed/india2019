all: poster

# PdfLaTeX compilation options
latexopt = -halt-on-error -file-line-error

poster: main.tex
	pdflatex $(latexopt) $< 
	pdflatex $(latexopt) $<
	pdflatex $(latexopt) $<
	mv main.pdf poster_fu.pdf

clean:
	  rm -f *.aux *.bbl *.blg *.log *.out *.snm *.toc *.nav

.PHONY: clean
