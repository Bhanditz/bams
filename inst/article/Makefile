DOC=HOCKING-model-selection-breakpoint-annotations

${DOC}.pdf: ${DOC}.tex refs.bib figure-1-smoothing.pdf figure-2-learning-curves.pdf algo.colors.tex figure-6-learned.png figure-roc2.pdf table-data-sets.tex figure-test-error-train-profiles.pdf figure-ann-sets.png figure-global-local.pdf figure-supplementary-18-anns.png table-unseen-detailed.tex figure-snp6-both.pdf table-training-error.tex global-original.tex 
	rm -f ${DOC}.aux ${DOC}.bbl
	pdflatex ${DOC}
	bibtex ${DOC}
	pdflatex ${DOC}
	pdflatex ${DOC}
segmentation.list.RData: segmentation.list.R
	R --no-save < $<
global-original.tex: global-original-detailed.R global.error.RData
	R --no-save < $<
table-training-error.tex: table-training-error.R profile.error.RData squares.R
	R --no-save < $<
figure-snp6.pdf: figure-snp6.R
	R --no-save < $<
figure-snp6-zoom.pdf: figure-snp6-zoom.R
	R --no-save < $<
figure-snp6-both.pdf: figure-snp6-both.R
	R --no-save < $<
table-unseen-detailed.tex: table-unseen-detailed.R error.unseen.RData timings.RData squares.R
	R --no-save < $<
timings.RData: timings.R segmentation.list.RData
	R --no-save < $<
figure-6-learned.png: figure-supplementary-18-anns.png
	cp figure-supplementary-207-anns.png figure-6-learned.png
figure-supplementary-18-anns.png: figure-supplementary.R signal.colors.R breakpoint.colors.R annotation.sets.RData
	R --no-save < $<
figure-global-local.pdf: figure-global-local.R leave.one.out.RData algo.colors.R
	R --no-save < $<
leave.one.out.RData: leave.one.out.R annotation.sets.RData fp.fn.list.RData
	R --no-save < $<
figure-ann-sets.png: figure-ann-sets.R breakpoint.colors.R
	R --no-save < $<
figure-test-error-train-profiles.pdf: figure-test-error-train-profiles.R error.unseen.RData algo.colors.R
	R --no-save < $<
error.unseen.RData: error.unseen.R profile.error.RData algos.in.tables.R
	R --no-save < $<
table-data-sets.tex: table-data-sets.R annotation.sets.RData
	R --no-save < $<
parameter.list.RData: parameter.list.R segmentation.list.RData
	R --no-save < $<
figure-roc2.pdf: figure-roc2.R global.error.RData algo.class.R algo.colors.R
	R --no-save < $<
global.error.RData: global.error.R profile.error.RData
	R --no-save < $<
profile.error.RData: profile.error.R fp.fn.list.RData annotation.sets.RData count.limits.R
	R --no-save < $<
fp.fn.list.RData: fp.fn.list.R segmentation.list.RData annotation.sets.RData count.limits.R
	R --no-save < $<
annotation.sets.RData: annotation.sets.R
	R --no-save < $<

figure-4-learned.png: figure-learned.R zzz.stats.RData signal.colors.R breakpoint.colors.R
	R --no-save < $<

table-annotation-chrom-counts.tex: table-annotation-chrom-counts.R
	R --no-save < $<

table-annotation-profile-counts.tex: table-annotation-profile-counts.R
	R --no-save < $<

table-error-on-unseen-profiles.tex: table-error-on-unseen-profiles.R zzz.stats.RData algos.in.tables.R
	R --no-save < $<

table-error-unseen-detailed.tex: table-error-unseen-detailed.R error.unseen.RData 
	R --no-save < $<
error.unseen.detailed.RData: error.unseen.detailed.R profile.error.RData algos.in.tables.R
	R --no-save < $<

table-generalization-error-global-models.tex: table-generalization-error-global-models.R zzz.stats.RData algos.in.tables.R
	R --no-save < $<

algo.colors.tex: algo.colors.R make.algo.colors.tex.R
	R --no-save < $<

figure-1-smoothing.pdf: figure-smoothing.R annotation.sets.RData
	R --no-save < $<

figure-2-learning-curves.pdf: figure-learning-curves.R breakpoint.colors.R profile.error.RData parameter.list.RData
	R --no-save < $<

figure-3-roc.pdf: figure-roc.R algo.colors.R zzz.stats.RData
	R --no-save < $<

figure-5-kinetics.pdf: figure-kinetics.R algo.colors.R zzz.stats.RData 
	R --no-save < $<

zzz.stats.RData: make.all.stats.R
	R --no-save < $<

clean:
	rm -rf table-*.tex *~ *.bbl *.blg *.log *.aux auto .Rhistory *.pdf *.out *.RData *.png