smoothdir <- file.path(Sys.getenv("HOME"),"smooth")
## TODO: criterion for done processing should be that we have as many
## algorithm subdirectories as the max
processed.cids <- dir(smoothdir)
counts <- sapply(processed.cids,function(cid){
  length(dir(file.path(smoothdir,cid)))
})
processed.cids <- names(counts)[counts==max(counts)]
library(bams)
data(neuroblastoma,package="bams")
all.cids <- levels(neuroblastoma$profiles$profile.id)
to.process <- all.cids[!all.cids %in% processed.cids]
algo <- "dnacopy.prune"
parameters <- eval(formals(smoothers[[algo]])[[2]])
profiles <- split(neuroblastoma$profiles,neuroblastoma$profiles$profile.id)
annotations <-
  split(neuroblastoma$annotations,neuroblastoma$annotations$profile.id)
for(cid in to.process){
  print(cid)
  breakpoint.labels <- annotations[[cid]]
  algodir <- file.path(smoothdir,cid,algo)
  dir.create(algodir)
  errors <- matrix(NA,nrow=length(parameters),ncol=nrow(breakpoint.labels))
  seconds <- Inf
  for(varname in c("breakpoint.labels","errors","parameters","seconds")){
    f <- gzfile(file.path(algodir,sprintf("%s.csv.gz",varname)),"w")
    write.table(get(varname),f,col.names=FALSE,row.names=FALSE,quote=FALSE,
                sep=",")
    close(f)
  }
}
