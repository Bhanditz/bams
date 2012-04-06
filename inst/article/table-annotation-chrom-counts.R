data(neuroblastoma,package="neuroblastoma")
library(plyr)
minmax <- ddply(neuroblastoma$annotations,.(chromosome),head,1)[,c("min","max")]
for(N in names(minmax)){
  minmax[[N]] <- sprintf("%.1f",minmax[[N]]/1e6)
}
colnames(minmax) <- c("min = $\\underline r_k$","max = $\\overline r_k$")
library(reshape2)
chrom.ann.counts <-
  dcast(neuroblastoma$annotations,chromosome~annotation,margins=TRUE)
colnames(chrom.ann.counts)[1] <- "chrom $c_k$"
library(xtable)
xt <- xtable(cbind(rbind(minmax,NA),chrom.ann.counts),align="rrrrrrr")
tex <- print(xt,include.rownames=FALSE,
             sanitize.text.function=identity,floating=FALSE)
FR <- list(c("hline","toprule"),
           c("hline","midrule"),
           c("hline","bottomrule"))
for(fr in FR){
  tex <- sub(fr[1],fr[2],tex)
}
cat(tex,file="table-annotation-chrom-counts.tex")
 
