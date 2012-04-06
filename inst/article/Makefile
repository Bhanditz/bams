DOC=HOCKING-model-selection-breakpoint-annotations

${DOC}.pdf: ${DOC}.tex refs.bib table-annotation-chrom-counts.tex table-generalization-error-global-models.tex figure-smoothing.pdf figure-learning-curves.pdf  figure-kinetics.pdf table-error-on-unseen-profiles.tex table-annotation-profile-counts.tex bioinfo.cls natbib.bst figure-roc.pdf algo.colors.tex
	rm -f ${DOC}.aux
	pdflatex ${DOC}
	bibtex ${DOC}
	pdflatex ${DOC}
	pdflatex ${DOC}

table-annotation-chrom-counts.tex: table-annotation-chrom-counts.R
	R --vanilla < table-annotation-chrom-counts.R

table-annotation-profile-counts.tex: table-annotation-profile-counts.R
	R --vanilla < table-annotation-profile-counts.R

table-error-on-unseen-profiles.tex: table-error-on-unseen-profiles.R zzz.stats.RData
	R --vanilla < table-error-on-unseen-profiles.R

table-generalization-error-global-models.tex: table-generalization-error-global-models.R zzz.stats.RData 
	R --vanilla < table-generalization-error-global-models.R

algo.colors.tex: algo.colors.R make.algo.colors.tex.R
	R --vanilla < make.algo.colors.tex.R

figure-roc.pdf: figure-roc.R algo.colors.R zzz.stats.RData
	R --vanilla < figure-roc.R

figure-smoothing.pdf: figure-smoothing.R 
	R --vanilla < figure-smoothing.R

figure-learning-curves.pdf: figure-learning-curves.R zzz.stats.RData
	R --vanilla < figure-learning-curves.R

figure-kinetics.pdf: figure-kinetics.R algo.colors.R zzz.stats.RData 
	R --vanilla < figure-kinetics.R

zzz.stats.RData: make.all.stats.R
	R --vanilla < make.all.stats.R

clean:
	rm -rf table-*.tex *~ *.bbl *.blg *.log *.aux auto .Rhistory figure-*.pdf all Rplots.pdf *.out
	touch zzz.stats.RData