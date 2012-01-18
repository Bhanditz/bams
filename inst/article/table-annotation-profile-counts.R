data(neuroblastoma,package="bams")
library(plyr)
counts <- ddply(neuroblastoma$annotations,.(profile.id),summarize,
      normal=sum(annotation=="normal"),
      breakpoints=sum(annotation=="breakpoint"))
## interesting summary table of counts per profile
counts.per.profile <- table(counts[,-1])
print(counts.per.profile)
library(xtable)
xt <- xtable(counts.per.profile,align="r|rrrrrrr")
tex <- print(xt,floating=FALSE)
tex <- sub("\\\\hline\n",
           paste("Normal&\\\\multicolumn{7}{c}{Breakpoint annotations}",
                 "\\\\\\\\",
                 "annotations"),tex)
cat(tex,file="table-annotation-profile-counts.tex")
