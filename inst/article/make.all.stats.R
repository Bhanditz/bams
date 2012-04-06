## the job of this script is to make zzz.stats.RData, which caches the
## results of the model smoothing in a nice R data structure
smoothdir <- file.path(Sys.getenv("HOME"),"smooth")
## TODO: maybe should actually check that all algo files are present
processed.cids <- dir(smoothdir)
nparams <- sapply(processed.cids,function(cid){
  algos <- dir(file.path(smoothdir,cid))
  sapply(algos,function(a){
    f <- file.path(smoothdir,cid,a,"parameters.csv.gz")
    if(file.exists(f)){
      length(scan(f,quiet=TRUE,what=1.2))
    }else 0
  })
},simplify=FALSE)
processed.algos <- unique(unlist(lapply(nparams,names)))
count.mat <- sapply(nparams,function(param.counts){
  param.counts[processed.algos]
})
## useful diagnostic == progress of the cluster
count.vecs <- apply(count.mat,1,table)
not.finished <- sapply(count.vecs,length)>1
print(count.vecs[not.finished])
## assume max is done...
done.mat <- apply(count.mat,1,function(x)x==max(x))
processed.cids <- rownames(done.mat)[apply(done.mat,1,all)]
data(neuroblastoma,package="neuroblastoma")
all.cids <- levels(neuroblastoma$profiles$profile.id)
to.process <- all.cids[!all.cids %in% processed.cids]
print(count.mat[not.finished,to.process])

## we have figured out which ones need processing, now do it.
profiles <- split(neuroblastoma$profiles,neuroblastoma$profiles$profile.id)
annotations <-
  split(neuroblastoma$annotations,neuroblastoma$annotations$profile.id)
library(bams)
for(cid in to.process){
  print(cid)
  done <- done.mat[cid,]
  algos.to.run <- names(done)[!done]
  pro <- profiles[[cid]]
  ann <- annotations[[cid]]
  run.smoothers(pro,ann,smoothers[algos.to.run])
}
algos <- dir(file.path(smoothdir,all.cids[1]))
all.stats <- list()
chrom.order <- as.character(c(1,2,3,4,11,17))
## each all.stats array is nparam x nprofiles x nann
for(a in algos){
  print(a)
  f <- file.path(smoothdir,processed.cids[1],a,"parameters.csv.gz")
  parameters <- scan(f,quiet=TRUE)
  param.names <- as.character(parameters)
  breakpoint.anns <-
    matrix(0,length(all.cids),length(chrom.order),
           dimnames=list(profile=all.cids,chromosome=chrom.order))
  normal.anns <- breakpoint.anns
  errors <-
    array(NA,
          list(length(param.names),length(all.cids),length(chrom.order)),
          list(param.names,all.cids,chrom.order))
  labels <- errors
  predictions <- errors
  false.positive <- errors
  false.negative <- errors
  for(cid in all.cids){
    errfile <- file.path(smoothdir,cid,a,"errors.csv.gz")
    e <- as.matrix(read.csv(errfile,header=FALSE))
    labfile <- file.path(smoothdir,cid,a,"breakpoint.labels.csv.gz")
    anns <-
      read.csv(labfile,header=FALSE,
               col.names=c("profile.id","chromosome","min","max","annotation"))
    ann.mat <- matrix(as.character(anns$annotation),
                      nrow=nrow(e),ncol=ncol(e),byrow=TRUE)
    colnames(e) <- anns$chromosome
    errors[,cid,colnames(e)] <- e
    for(chr in colnames(e)){
      normal.anns[cid,chr] <- if(is.na(e[1,chr]))0 else
        nrow(subset(anns,chromosome==chr & annotation=="normal"))
      breakpoint.anns[cid,chr] <- if(is.na(e[1,chr]))0 else
        nrow(subset(anns,chromosome==chr & annotation=="breakpoint"))
    }
    ## Careful: is.na(NA & TRUE) but !is.na(NA & FALSE)
    false.negative[,cid,colnames(e)] <-
      ifelse(is.na(e),NA,ifelse(e & ann.mat=="breakpoint",1L,0L))
    false.positive[,cid,colnames(e)] <-
      ifelse(is.na(e),NA,ifelse(e & ann.mat=="normal",1L,0L))
  }
  readsecs <- function(cid){
    secfile <- file.path(smoothdir,cid,a,"seconds.csv.gz")
    scan(secfile,quiet=TRUE)
  }
  seconds <- sapply(all.cids,readsecs)
  all.stats[[a]] <- list(errors=errors,
                         false.positive=false.positive,
                         false.negative=false.negative,
                         parameters=parameters,
                         seconds=seconds,
                         normal.anns=normal.anns,
                         breakpoint.anns=breakpoint.anns)
}
save(all.stats,file="zzz.stats.RData")
