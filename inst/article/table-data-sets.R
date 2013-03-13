works_with_R("2.15.2",xtable="1.7.0")
load("annotation.sets.RData")
sets <- names(annotation.sets)
anns <- c("0breakpoints","1breakpoint",">0breakpoints")
protocols <- c(original="Systematic",detailed="Any")
info <- do.call(cbind,lapply(sets,function(set.name){
  set <- annotation.sets[[set.name]]
  counts <- sapply(anns,function(a){
    sum(set$ann==a)
  })
  c(protocol=protocols[[set.name]],
    "annotated profiles"=length(unique(set[,"profile.id"])),
    "annotated chromosomes"=nrow(unique(set[,c("profile.id","chromosome")])),
    annotations=nrow(set),
    counts)
}))
colnames(info) <- sets
xt <- xtable(info,align="rrr")
rownames(xt) <- sub(">","$>$",rownames(xt))
print(xt,floating=FALSE,sanitize.text.function=identity,
      file="table-data-sets.tex")
