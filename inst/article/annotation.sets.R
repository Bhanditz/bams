data(neuroblastoma,package="neuroblastoma")
data(neuroblastomaDetailed,package="bams")

annotation.sets <- list(original=neuroblastoma$annotations,
                        detailed=neuroblastomaDetailed)

## standardize annotation levels.
standard <- c(breakpoint=">0breakpoints",
              normal="0breakpoints")
for(set.name in names(annotation.sets)){
  nonstandard <- annotation.sets[[set.name]]$ann%in%names(standard)
  annotation.sets[[set.name]]$annotation <- 
    as.character(annotation.sets[[set.name]]$annotation)
  annotation.sets[[set.name]][nonstandard,"annotation"] <-
    standard[annotation.sets[[set.name]][nonstandard,"annotation"]]
  stopifnot(all(!annotation.sets[[set.name]]$ann%in%names(standard)))
}

save(annotation.sets,file="annotation.sets.RData")
