all:
	pdflatex draft.tex
	bibtex draft
	pdflatex draft.tex
	pdflatex draft.tex

clean:
	@rm -f draft.pdf draft.aux draft.blg draft.bbl
	@rm -f draft.log
