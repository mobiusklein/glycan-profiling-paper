.SECONDEXPANSION:


full: clean
	pdflatex draft.tex
	bibtex draft
	pdflatex supplement.tex
	bibtex supplement
	pdflatex draft.tex
	pdflatex supplement.tex
	pdflatex draft.tex
	pdflatex supplement.tex

draft: $(FIGURES) supplement.pdf
	pdflatex draft.tex
	bibtex draft
	pdflatex draft.tex
	pdflatex draft.tex

layout:
	pdflatex draft.tex

sup:
	pdflatex supplement.tex
	bibtex supplement
	pdflatex supplement.tex
	pdflatex supplement.tex


clean:
	@rm -f draft.pdf draft.aux draft.blg draft.bbl draft.log draft.out
	@rm -f supplement.pdf supplement.aux supplement.blg supplement.bbl supplement.log supplement.out

FIGURES =
