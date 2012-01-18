data(neuroblastoma,package="bams")
clin.id <- commandArgs(trailingOnly=TRUE)
one <- subset(neuroblastoma$profiles,profile.id==clin.id)
these.labels <- subset(neuroblastoma$annotations,profile.id==clin.id)
print(these.labels)
library(bams)
run.smoothers(one,these.labels,article.smoothers)

