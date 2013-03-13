works_with_R("2.15.2",xtable="1.7.0",plyr="1.8")

load("error.unseen.RData")
load("timings.RData")
source("squares.R")

detailed <- subset(error.unseen,
                   train.size==10 &
                   train.name=="detailed" &
                   test.name=="detailed")

results <- ddply(detailed,.(algo,type),summarize,
                 mean=mean(error),sd=sd(error))

rlist <- split(results,results$type)

mean.sd <- do.call(cbind,lapply(rlist,"[",c("mean","sd")))
ord <- order(mean.sd$error.mean)
## do not re-order yet!!
mean.sd <- sapply(mean.sd,function(x)sprintf("%.1f",x*100))
## add colored squares
algos <- as.character(rlist[[1]]$algo)
square.inches <- "0.08"
errsd <-
  data.frame(algorithm=squares(algos),
             mean.sd,
             Timings=sprintf("%.2f",timings[algos],2))

## nrow(all.est[[1]]) ## number of folds
## ## TODO: Correct FP/FN. DONT divide by total number of annotations!
## errtab <- t(sapply(all.est,colMeans))
## errtab <- errtab[order(errtab[,1]),]
## sdtab <- t(sapply(all.est,function(x)apply(x,2,sd)))[rownames(errtab),]
## errsd <- do.call(cbind,lapply(colnames(errtab),function(i){
##   x <- cbind(errtab[,i],sdtab[,i])
##   colnames(x) <- c(i,"sd")
##   x
## }))
## colnames(errsd)[c(3,5)] <- c("FP","FN")


xt <- xtable(errsd[ord,],align="rrrrrrrrr")
colnames(xt)[c(3,5,7)] <- "sd"
colnames(xt) <- sub(".mean","",colnames(xt))
tex <- print(xt,floating=FALSE,include.rownames=FALSE,
             sanitize.text.function=identity)
cat(tex,file="table-unseen-detailed.tex")
