works_with_R("2.15.2",bams="1.6")

load("profile.error.RData")
max.train.size <- 30
pids <- rownames(profile.error[[1]][[1]][[1]])
algos <- c("glad.lambdabreak","cghseg.k","flsa.norm","dnacopy.sd")
annotators <- data.frame()
get.err <- function(L,pids){
  with(L,colSums((fp+fn)[pids,,drop=FALSE],na.rm=TRUE))
}
for(algo in algos){
  for(train.size in 1:max.train.size){
    cat(sprintf("%3d / %3d train %20s\n",
                train.size,max.train.size,algo))
    nfolds <- floor(length(pids)/train.size)
    fold <- sample(rep(1:nfolds,l=length(pids)))
    ## first split on folds.
    for(f in 1:nfolds){
      is.test <- fold != f
      is.train <- !is.test
      train.pids <- pids[is.train]
      test.pids <- pids[is.test]
      ## then for every fold do a 2 x 2 training/testing on one data
      ## set or the other.
      for(train.name in names(profile.error)){
        train.err <- get.err(profile.error[[train.name]][[algo]],train.pids)
        picked <- pick.best.index(train.err)
        for(test.name in names(profile.error)){
          errors <-
            get.err(profile.error[[test.name]][[algo]],test.pids)[picked]
          anns <- sum(profile.error[[test.name]][[algo]]$anns[test.pids,])
          annotators <- rbind(annotators,{
            data.frame(algo,train.size,train.name,test.name,error=errors/anns)
          })
        }
      }
    }
  }
}

save(annotators, file="annotators.RData")
