load("segmentation.list.RData")
load("annotation.sets.RData")
source("count.limits.R")

algos <- names(segmentation.list[[1]])

fp.fn.list <- list()
for(set.name in names(annotation.sets)){
  set <- annotation.sets[[set.name]]
  for(algo in algos){
    params <- segmentation.list[[1]][[algo]]$parameters
    result <- list(fp=matrix(NA,nrow(set),length(params)),
                   fn=matrix(NA,nrow(set),length(params)),
                   tp=matrix(NA,nrow(set),length(params)))
    for(i in 1:nrow(set)){
      cat(sprintf("%20s %10s %5d / %5d\n",algo,set.name,i,nrow(set)))
      a <- set[i,]
      ## counts of breakpoints in this window for every parameter of
      ## this algorithm.
      pid <- as.character(a$profile.id)
      chr <- as.character(a$chromosome)
      seg.info <- segmentation.list[[pid]][[algo]]
      if(is.null(seg.info)){
        warning(sprintf("no segmentation for %s %s, NAs generated",
                        pid, algo))
      }else{
        counts <- sapply(seg.info$breakpoints[[chr]],function(positions){
          sum(a$min < positions & positions < a$max)
        })
        ann <- as.character(a$annotation)
        lim <- count.limits[[ann]]
        if(is.null(lim)){
          warning(sprintf("undefined cost for %s, NAs generated",ann))
        }else{
          result$fp[i,] <- counts > lim[2]
          result$fn[i,] <- counts < lim[1]
          result$tp[i,] <- counts >= lim[1] & lim[1] > 0
        }
      }
    }
    fp.fn.list[[set.name]][[algo]] <- result
  }
}
        
save(fp.fn.list, file="fp.fn.list.RData")
