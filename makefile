# Name of manuscript
manuscript = main

.PHONY: $(manuscript).pdf all clean

all: $(manuscript).pdf

#=================================================================
# Generate PDF of manuscript using PdfLaTeX
#=================================================================

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interaction=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

main.pdf: main.tex
	latexmk -pdf -pdflatex="pdflatex -interaction=nonstopmode" -use-make $(manuscript).tex

clean:
	latexmk -CA
