works_with_R("2.15.0",gada="1.0",changepoint="1.0.4",
             DNAcopy="1.18.0",cghseg="0.1",flsa="1.3",
             cghFLasso="0.2.1",GLAD="2.0.0",bams="1.6",
             neuroblastoma="1.0")

data(neuroblastoma)
clin.id <- commandArgs(trailingOnly=TRUE)
profile <- subset(neuroblastoma$profiles,profile.id==clin.id)
print(head(profile))
library(bams)
#update.names <- grep("iir",names(smoothers),value=TRUE,invert=TRUE)
#update.smoothers <- article.smoothers
update.smoothers <- smoothers["cghseg.k"]
seg.profile(profile, smooth.funs=update.smoothers)

