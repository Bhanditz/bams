data(neuroblastoma,package="neuroblastoma")
clin.id <- commandArgs(trailingOnly=TRUE)
one <- subset(neuroblastoma$profiles,profile.id==clin.id)
these.labels <- subset(neuroblastoma$annotations,profile.id==clin.id)
print(these.labels)
library(bams)
update.names <- c("gada",grep("cghseg",names(smoothers),value=TRUE))
update.smoothers <- smoothers[update.names]
run.smoothers(one, these.labels, update.smoothers)

