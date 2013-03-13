load("segmentation.list.RData")

anames <- names(segmentation.list[[1]])

sec.mat <- t(sapply(segmentation.list,function(algo.list){
  sapply(algo.list,function(L){
    L$seconds
  })[anames]
})) ## pid x algo
sec.mat[is.na(sec.mat)] <- Inf

timings <- sort(apply(sec.mat, 2, median))

save(timings, file="timings.RData")
