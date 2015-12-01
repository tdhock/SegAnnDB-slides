SLIDES=HOCKING-interactive-genomic-segmentation
${SLIDES}.pdf: ${SLIDES}.tex figure-profiles.png figure-annotations.png figure-segannot-1.pdf
	rm -f *.aux *.bbl
	pdflatex ${SLIDES}
	pdflatex ${SLIDES}
figure-profiles.png: figure-profiles.R clinical-limited.csv
	R --no-save < $<
signal.list.RData: signal.list.R 
	R --no-save < $<
segmentation.list.RData: segmentation.list.R signal.list.RData
	R --no-save < $<
signal.features.RData: signal.features.R signal.list.RData
	R --no-save < $<
exact.breakpoints.RData: exact.breakpoints.R segmentation.list.RData
	R --no-save < $<
figure-annotations.png: figures-demo.R signal.list.RData demo.csv segmentation.list.RData exact.breakpoints.RData signal.features.RData
	R --no-save < $<
demo:
	rm -f demo.csv
	make
demo.csv: profiles.csv
	python annotate_breakpoints.py profiles.csv demo.csv
profiles.csv: make-profiles.R clinical-limited.csv
	R --no-save < $<
figure-segannot-1.pdf: figure-segannot.R
	R --no-save < $<
.PHONY: clean
clean: 
	rm *.RData