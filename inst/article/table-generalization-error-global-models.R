estimate.test.error <- function
### Do leave-one-out cross-validation on chromosome arms.
(stats
### Named list with arrays errors, false.positive, false.negative,
### each of dim nparam x nprof x nfolds.
 ){
  stats <- stats[c("errors","false.positive","false.negative")]
  nparam <- dim(stats$errors)[1]
  nprof <- dim(stats$errors)[2]
  nfolds <- dim(stats$errors)[3]
  stat.names <- names(stats)
  local.loo <- matrix(NA,length(stats),nfolds)
  hybrid.loo <- matrix(NA,length(stats),nfolds)
  global.loo <- matrix(NA,length(stats),nfolds)
  rownames(global.loo) <- stat.names
  rownames(hybrid.loo) <- stat.names
  rownames(local.loo) <- stat.names
  train.err.mat <-
    matrix(NA,3,nfolds,dimnames=list(method=c("global","hybrid","local"),fold=NULL))
  for(fold in 1:nfolds){

    train.err <- rep(NA,nparam) ## global model
    for(j in 1:nparam){
      train.err[j] <- sum(stats$errors[j,,-fold],na.rm=TRUE)
    }
    ## save for hybrid approach
    global.train.err <-
      data.frame(train.err,param=1:length(train.err))
    global.picked <- pick.best.index(train.err)
    for(sn in stat.names){
      global.loo[sn,fold] <-
        mean(stats[[sn]][global.picked,,fold],na.rm=TRUE)
    }
    train.err.mat["global",fold] <- train.err[global.picked] ## for comparing train err

    ind.stats <- matrix(NA,nprof,length(stats))
    colnames(ind.stats) <- stat.names
    hybrid.train.errors <- rep(NA,nprof)
    for(i in 1:nprof){ ## hybrid models
      train.err <-
        apply(stats$errors[,i,-fold,drop=FALSE],1,sum,na.rm=TRUE)
      is.min <- train.err == min(train.err)
      global.subset <- global.train.err[is.min,]
      hybrid.picked <-
        with(global.subset,param[which.min(train.err)])
      hybrid.train.errors[i] <- train.err[hybrid.picked]
      for(sn in stat.names){ ## store test err for picked model
        ind.stats[i,sn] <- stats[[sn]][hybrid.picked,i,fold]
      }
    }
    hybrid.loo[,fold] <- colMeans(ind.stats,na.rm=TRUE)
    train.err.mat["hybrid",fold] <- sum(hybrid.train.errors)
    
    ind.stats <- matrix(NA,nprof,length(stats))
    colnames(ind.stats) <- stat.names
    local.train.errors <- rep(NA,nprof)
    for(i in 1:nprof){ ## local models
      train.err <-
        apply(stats$errors[,i,-fold,drop=FALSE],1,sum,na.rm=TRUE)
      local.picked <- pick.best.index(train.err)
      for(sn in stat.names){ ## store test err for picked model
        ind.stats[i,sn] <- stats[[sn]][local.picked,i,fold]
      }
      local.train.errors[i] <- train.err[local.picked]
    }
    local.loo[,fold] <- colMeans(ind.stats,na.rm=TRUE)
    train.err.mat["local",fold] <- sum(local.train.errors)
    
  }
  list(local=local.loo,
       hybrid=hybrid.loo,
       global=global.loo,
       train.err.mat=train.err.mat)
### Named list with elements local, hybrid, global, each a 3 x nfolds
### matrix.
}


load("zzz.stats.RData")
library(bams)
source("algos.in.tables.R")
all.stats <- all.stats[algos.in.tables]
cghseg.k.result <- estimate.test.error(all.stats$cghseg.k)
cghseg.k.result$train
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
library(xtable)
## html code:
## header <- readLines("header.html")
## footer <- readLines("footer.html")
## html.tab <- fpfn
## html.xt <- xtable(html.tab,align="rrrrrrrr")
## html <- print(html.xt,type="html",html.table.attributes="")
## find <- "<TABLE >"
## firstlines <- paste(c(find,
##                       "<tr>",
##                       "<td></td>",
##                       "<th align=center colspan=3>global</th>",
##                       "<th align=center colspan=3>local</th>",
##                       "<th>Timings</th>",
##                       "</tr>"),collapse="")
## edited <- gsub(find,firstlines,html,fixed=TRUE)
## with.header <- c(header,edited,footer)
## cat(with.header,file="~/public_html/neuroblastoma/accuracy.html")

## tex code:
counts.percents <- fpfn
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
