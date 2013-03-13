load("fp.fn.list.RData")
load("annotation.sets.RData")
source("count.limits.R")

## list[[ann.set]][[algo]] is a list of 2 matrices: fp and fn, each of
## which is n.annotations x n.parameters. The job of this script is to
## convert that into local error functions for each profile:
## list[[ann.set]][[algo]][[fp/fn/anns]][profile.id,parameter].

pids <- levels(annotation.sets[[1]]$profile.id)

pid <- "520"
algo <- "dnacopy.prune"

profile.error <- list()
for(set.name in names(annotation.sets)){
  anns <- annotation.sets[[set.name]]
  err <- fp.fn.list[[set.name]]
  anns.by.pid <- split(anns,anns$profile.id)
  for(algo in names(err)){
    cat(sprintf("%10s %20s\n",set.name,algo))
    result <- list(anns=matrix(NA,length(pids),length(count.limits),
                   dimnames=list(pids,names(count.limits))))
    for(stat in names(err[[algo]])){
      result[[stat]] <- matrix(NA,length(pids),ncol(err[[algo]]$fp))
      rownames(result[[stat]]) <- pids
    }
    df.list <- lapply(err[[algo]],function(m){
      split(data.frame(m),anns$profile.id)
    })
    for(pid in pids){
      for(var.name in names(err[[algo]])){
        result[[var.name]][pid,] <-
          sapply(df.list[[var.name]][[pid]],sum,na.rm=TRUE)
      }
      ann.f <- factor(anns.by.pid[[pid]]$annotation,names(count.limits))
      result$anns[pid,] <-
        as.integer(table(ann.f[!is.na(df.list$fp[[pid]][[1]])]))
    }
    profile.error[[set.name]][[algo]] <- result
  }
}

save(profile.error, file="profile.error.RData")
