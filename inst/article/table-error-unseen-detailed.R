works_with_R("2.15.2",xtable="1.7.0",plyr="1.8")

load("error.unseen.RData")

detailed <- subset(error.unseen,
                   train.size==10 &
                   train.name=="detailed" &
                   test.name=="detailed")

results <- ddply(detailed,.(algo,type),summarize,
                 mean=mean(error),sd=sd(error))

rlist <- split(results,results$type)

mean.sd <- do.call(cbind,lapply(rlist,"[",c("mean","sd")))

nrow(all.est[[1]]) ## number of folds
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

## add colored squares
square.inches <- "0.08"
rownames(errsd) <-
  sprintf("%s \\textcolor{%s.color}{\\rule{%sin}{%sin}}",
          rownames(errsd),
          rownames(errsd),
          square.inches,
          square.inches)

xt <- xtable(data.frame(model=rownames(errsd),errsd*100),
             digits=1,align="rrrrrrrr")
colnames(xt)[c(3,5,7)] <- "sd"
tex <- print(xt,floating=FALSE,include.rownames=FALSE,
             sanitize.text.function=identity)
##tex <- sub("cghseg.mBIC","\\\\hline cghseg.mBIC",tex)
FR <- list(c("hline","toprule"),
           c("hline","midrule"),
           c("hline","midrule"),
           c("hline","bottomrule"))
## for(fr in FR){
##   tex <- sub(fr[1],fr[2],tex)
## }
cat(tex,file="table-error-on-unseen-profiles.tex")
