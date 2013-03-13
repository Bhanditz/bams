works_with_R("2.15.2",bams="1.6",neuroblastoma="1.0")

## the job of this script is to read the data saved by the smoothing
## algorithms to the ~/seg directory.
segdir <- file.path(Sys.getenv("HOME"),"seg")
data(neuroblastoma)
processed.cids <- as.character(unique(neuroblastoma$annotations$pro))

## see if pelt.n is processed new-style.
first <- sapply(processed.cids,function(cid){
  f <- file.path(segdir,cid,"pelt.n","parameters.csv.gz")
  if(file.exists(f)){
    scan(f,nmax=1,quiet=TRUE,what="char")
  }else NA
})
not <- names(first)[grepl("n",first)]

## TODO: maybe should actually check that all algo files are present
nparams <- sapply(processed.cids,function(cid){
  algos <- dir(file.path(segdir,cid))
  sapply(algos,function(a){
    f <- file.path(segdir,cid,a,"parameters.csv.gz")
    if(file.exists(f)){
      length(scan(f,quiet=TRUE,what="char"))
    }else 0
  })
},simplify=FALSE)
processed.algos <- unique(unlist(lapply(nparams,names)))
count.mat <- sapply(nparams,function(param.counts){
  param.counts[processed.algos]
})
count.mat[is.na(count.mat)] <- 0
## useful diagnostic == progress of the cluster
count.vecs <- apply(count.mat,1,table)
not.finished <- sapply(count.vecs,length)>1
print(count.vecs[not.finished])
## assume max is done...
done.mat <- apply(count.mat,1,function(x)x==max(x))
processed.cids <- rownames(done.mat)[apply(done.mat,1,all)]
all.cids <- levels(neuroblastoma$profiles$profile.id)
to.process <- all.cids[!all.cids %in% processed.cids]
print(count.mat[not.finished,to.process,drop=FALSE])

## we have figured out which ones need processing, now do it.
profiles <- split(neuroblastoma$profiles,neuroblastoma$profiles$profile.id)
for(cid in to.process){
  done <- done.mat[cid,]
  algos.to.run <- names(done)[!done]
  #algos.to.run <- grep("dnacopy.prune",algos.to.run,invert=TRUE,value=TRUE)
  pro <- profiles[[cid]]
  cat(sprintf("profile %s, %d probes\n",pro$pro[1],nrow(pro)))
  seg.profile(pro, smoothers[algos.to.run])
}

read.algo <- function(algo.dir){
  results <- list()
  getf <- function(x)sprintf("%s/%s.csv.gz",algo.dir,x)
  for(base in c("parameters","seconds")){
    results[[base]] <- scan(getf(base),quiet=TRUE)
  }
  break.df <-
    read.csv(getf("breakpoints"),header=FALSE,
             col.names=c("parameter","chromosome","base"))
  break.df$chromosome <- factor(break.df$chromosome,c(1:22,"X"))
  break.df$parameter <- factor(break.df$parameter,results$parameters)
  results$breakpoints <- lapply(split(break.df,break.df$chrom),function(d){
    df.list <- split(d,d$parameter)
    lapply(df.list, "[[", "base")
  })
  results
}
  
read.profile <- function(pid.dir){
  results <- list()
  algos <- dir(pid.dir)
  for(algo in algos){
    results[[algo]] <- read.algo(file.path(pid.dir, algo))
  }
  results
}

pid.dirs <- Sys.glob(file.path(segdir,"*"))
segmentation.list <- list()
for(i in seq_along(pid.dirs)){
  pdir <- pid.dirs[[i]]
  cat(sprintf("%4d / %4d %20s\n",i, length(pid.dirs), pdir))
  segmentation.list[[pid]] <- read.profile(pdir)
}

save(segmentation.list, file="segmentation.list.RData")
