.SECONDEXPANSION:


all: $(FIGURES) supplement.pdf
	# pdflatex -draftmode supplement.tex
	pdflatex -draftmode draft.tex
	bibtex draft
	pdflatex -draftmode draft.tex
	pdflatex draft.tex


supplement:
	pdflatex -draftmode supplement.tex
	bibtex supplement
	pdflatex -draftmode supplement.tex
	pdflatex supplement.tex


clean:
	@rm -f draft.pdf draft.aux draft.blg draft.bbl
	@rm -f draft.log

FIGURES =

figure/igg_chromatograms.pdf:
	inkscape -z -f "figure/igg_chromatograms.svg" -A "figure/igg_chromatograms.pdf"
FIGURES += figure/igg_chromatograms.pdf

figure/igg_abundances.pdf:
	inkscape -z -f "figure/igg_abundances.svg" -A "figure/igg_abundances.pdf"	
FIGURES += figure/igg_abundances.pdf

figure/phil_82_chromatograms.pdf:
	inkscape -z -f "figure/phil_82_chromatograms.svg" -A "figure/phil_82_chromatograms.pdf"
FIGURES += figure/phil_82_chromatograms.pdf

figure/phil_82_abundances.pdf:
	inkscape -z -f "figure/phil_82_abundances.svg" -A "figure/phil_82_abundances.pdf"
FIGURES += figure/phil_82_abundances.pdf

figure/agp_chromatograms.pdf:
	inkscape -z -f "figure/agp_chromatograms.svg" -A "figure/agp_chromatograms.pdf"
FIGURES += figure/agp_chromatograms.pdf

figure/agp_abundances.pdf:
	inkscape -z -f "figure/agp_abundances.svg" -A "figure/agp_abundances.pdf"	
FIGURES += figure/agp_abundances.pdf

figure/phil_bs_chromatograms.pdf:
	inkscape -z -f "figure/phil_bs_chromatograms.svg" -A "figure/phil_bs_chromatograms.pdf"
FIGURES += figure/phil_bs_chromatograms.pdf

figure/phil_bs_abundances.pdf:
	inkscape -z -f "figure/phil_bs_abundances.svg" -A "figure/phil_bs_abundances.pdf"	
FIGURES += figure/phil_bs_abundances.pdf

figure/dp_phil_82_chromatograms.pdf:
	inkscape -z -f "figure/dp_phil_82_chromatograms.svg" -A "figure/dp_phil_82_chromatograms.pdf"
FIGURES += figure/dp_phil_82_chromatograms.pdf

figure/dp_phil_82_abundances.pdf:
	inkscape -z -f "figure/dp_phil_82_abundances.svg" -A "figure/dp_phil_82_abundances.pdf"	
FIGURES += figure/dp_phil_82_abundances.pdf

figure/dp_agp_chromatograms.pdf:
	inkscape -z -f "figure/dp_agp_chromatograms.svg" -A "figure/dp_agp_chromatograms.pdf"
FIGURES += figure/dp_agp_chromatograms.pdf

figure/dp_agp_abundances.pdf:
	inkscape -z -f "figure/dp_agp_abundances.svg" -A "figure/dp_agp_abundances.pdf"	
FIGURES += figure/dp_agp_abundances.pdf

figure/rp_agp_chromatograms.pdf:
	inkscape -z -f "figure/rp_agp_chromatograms.svg" -A "figure/rp_agp_chromatograms.pdf"
FIGURES += figure/rp_agp_chromatograms.pdf

figure/rp_agp_abundances.pdf:
	inkscape -z -f "figure/rp_agp_abundances.svg" -A "figure/rp_agp_abundances.pdf"	
FIGURES += figure/rp_agp_abundances.pdf

figure/rp_human_serum_chromatograms.pdf:
	inkscape -z -f "figure/rp_human_serum_chromatograms.svg" -A "figure/rp_human_serum_chromatograms.pdf"
FIGURES += figure/rp_human_serum_chromatograms.pdf

figure/rp_human_serum_abundances.pdf:
	inkscape -z -f "figure/rp_human_serum_abundances.svg" -A "figure/rp_human_serum_abundances.pdf"	
FIGURES += figure/rp_human_serum_abundances.pdf


figure/process_flowchart.pdf:
	inkscape -z -f "figure/process_flowchart.svg" -A "figure/process_flowchart.pdf"	
FIGURES += figure/process_flowchart.pdf


figures: $(FIGURES)
	@echo $(FIGURES)

clean-figures:
	rm $(FIGURES)
