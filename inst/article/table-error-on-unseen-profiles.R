works_with_R("2.15.2",xtable="1.7.0")

unseen.profile.error <- function
### Do n/t-fold cross-validation to estimate the error of global
### models with a small training set.
(stats,
### list with arrays for errors, false.positive, and false.negative
 prof.per.train
### t = approximate number of annotated profiles per training set.
 ){
  stat.names <- c("errors","false.positive","false.negative")
  nparam <- dim(stats$errors)[1]
  nprof <- dim(stats$errors)[2]
  narms <- dim(stats$errors)[3]
  nfolds <- floor(nprof/prof.per.train)
  test.error <- matrix(NA,nfolds,length(stat.names))
  colnames(test.error) <- stat.names
  fold <- sample(rep(1:nfolds,l=nprof))
  for(f in 1:nfolds){
    test.fold <- fold != f
    train.fold <- !test.fold
    train.err <- rep(NA,nparam) ## global model
    for(j in 1:nparam){
      train.err[j] <- mean(stats$errors[j,train.fold,],na.rm=TRUE)
    }
    global.picked <- pick.best.index(train.err)
    for(sn in stat.names){
      test.error[f,sn] <-
        mean(stats[[sn]][global.picked,test.fold,],na.rm=TRUE)
    }
    num.normal <- sum(stats$normal.anns[test.fold,])
    num.breakpoint <- sum(stats$breakpoint.anns[test.fold,])
    num.total <- num.normal+num.breakpoint
    FP <- test.error[f,"false.positive"]*num.total
    if(round(FP)-FP>1e-5){
      stop("FP not integral!")
    }
    test.error[f,"false.positive"] <- FP/num.normal
    test.error[f,"false.negative"] <-
      test.error[f,"false.negative"]*num.total/num.breakpoint
  }
  test.error
### matrix of estimated test errors, nfolds x 3
}

load("zzz.stats.RData")
all.est <- list()
source("algos.in.tables.R")
for(algo in algos.in.tables){
  set.seed(2)
  all.est[[algo]] <- unseen.profile.error(all.stats[[algo]],10)
}
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
