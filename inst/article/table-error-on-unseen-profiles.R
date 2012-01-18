library(bams)
load("zzz.stats.RData")
all.est <- list()
for(algo in names(all.stats)){
  set.seed(2)
  all.est[[algo]] <- unseen.profile.error(all.stats[[algo]],10)
}
## TODO: Correct FP/FN. DONT divide by total number of annotations!
errtab <- t(sapply(all.est,colMeans))
errtab <- errtab[order(errtab[,1]),]
sdtab <- t(sapply(all.est,function(x)apply(x,2,sd)))[rownames(errtab),]
errsd <- do.call(cbind,lapply(colnames(errtab),function(i){
  x <- cbind(errtab[,i],sdtab[,i])
  colnames(x) <- c(i,"sd")
  x
}))
colnames(errsd)[c(3,5)] <- c("FP","FN")
library(xtable)
no.tuning <- rownames(errsd)%in%
  c("glad.default","dnacopy.default","cghseg.mBIC","cghFLasso")
rownames(errsd)[no.tuning] <- paste("*",rownames(errsd)[no.tuning],sep="")
xt <- xtable(data.frame(model=rownames(errsd),errsd*100),
             digits=1,align="rrrrrrrr")
colnames(xt)[c(3,5,7)] <- "sd"
tex <- print(xt,floating=FALSE,include.rownames=FALSE)
##tex <- sub("cghseg.mBIC","\\\\hline cghseg.mBIC",tex)
FR <- list(c("hline","toprule"),
           c("hline","midrule"),
           c("hline","midrule"),
           c("hline","bottomrule"))
## for(fr in FR){
##   tex <- sub(fr[1],fr[2],tex)
## }
cat(tex,file="table-error-on-unseen-profiles.tex")
