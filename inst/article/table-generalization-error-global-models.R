load("zzz.stats.RData")
library(bams)
## OK: MBIC is the same for global and local.
if("cghseg.mBIC" %in% names(all.stats)){
  print(with(estimate.test.error(all.stats$cghseg.mBIC),local-global))
}
results <- lapply(all.stats,estimate.test.error)
## compare local, hybrid, and global error rates
print(sapply(results,function(stat)sapply(stat,function(x)mean(x[1,]))))
library(plyr)
arm.generalization.error <- laply(results,"[[","global")
names(dimnames(arm.generalization.error)) <- c("algorithm","statistic","chromosome")
dimnames(arm.generalization.error)[[1]] <- names(all.stats)


## average over all arms
global.local.stats <- laply(results,function(L)do.call(rbind,lapply(L,rowMeans)))
names(dimnames(global.local.stats)) <- c("algorithm","method","statistic")
dimnames(global.local.stats)[[1]] <- names(all.stats)
dimnames(global.local.stats)[[3]] <- c("errors","FP","FN")

## now figure out a good ordering
arm.error <- arm.generalization.error[,1,]
ord <- rownames(arm.error)[order(apply(arm.error,1,mean))]

## old version (divide by total number of examples)
fpfn <- cbind(global.local.stats[,"global",],
              global.local.stats[,"local",])
## new version (divide by number of positive/neg examples)
num.normal <- do.call(rbind,lapply(all.stats,function(stat){
  colSums(stat$normal.anns)
}))
num.breakpoint <- do.call(rbind,lapply(all.stats,function(stat){
  colSums(stat$breakpoint.anns)
}))
num.total <- num.normal+num.breakpoint
fp.count <- arm.generalization.error[,"false.positive",]*num.total
stopifnot(all(abs(round(fp.count)-fp.count)<1e-6))
fp <- fp.count/num.normal
fn <- arm.generalization.error[,"false.negative",]*num.total/num.breakpoint
## do for each local and global
fpfn.list <- list()
for(training.method in c("global","local")){
  errors <- rep(NA,length(results))
  names(errors) <- names(results)
  FP <- rep(NA,length(results))
  names(FP) <- names(results)
  FN <- rep(NA,length(results))
  names(FN) <- names(results)
  for(algorithm in names(results)){
    m <- results[[algorithm]][[training.method]]
    FP[algorithm] <-
      mean(m["false.positive",]*num.total[algorithm,]/num.normal[algorithm,])
    FN[algorithm] <-
     mean(m["false.negative",]*num.total[algorithm,]/num.breakpoint[algorithm,])
    errors[algorithm] <- mean(m["errors",])
  }
  fpfn.list[[training.method]] <- cbind(errors,FP,FN)
}
fpfn <- with(fpfn.list,cbind(global,local))[ord,]*100
fpfn[1:length(fpfn)] <- sprintf("%3.1f",fpfn)


## add timings
seconds <- lapply(all.stats,"[[","seconds")
fpfn <- cbind(fpfn,seconds=sprintf("%.2f",sapply(seconds,median)[ord]))

counts.percents <- fpfn
library(xtable)
square.inches <- "0.08"
rownames(counts.percents) <-
  sprintf("%s \\textcolor{%s.color}{\\rule{%sin}{%sin}}",
          rownames(counts.percents),
          rownames(counts.percents),
          square.inches,
          square.inches)
        
xt <- xtable(counts.percents,align="r|rrr|rrr|r")
tex <- print(xt,floating=FALSE,sanitize.rownames.function=identity,
             include.rownames=TRUE,size="small")
REP <- paste("\\hline &",
             "\\\\multicolumn{3}{c|}{Global} &",
             "\\\\multicolumn{3}{c|}{Local}  &",
             "Timings ",
             "\\\\\\\\")
tex <- sub("\\hline",REP,tex)
##tex <- sub("cghseg.mBIC","\\\\hline cghseg.mBIC",tex)
tex <- sub("\\hline(\nannotations.*?\n)","\\1\\\\hline",tex)
cat(tex,file="table-generalization-error-global-models.tex")
