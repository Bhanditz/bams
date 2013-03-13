works_with_R("2.15.2",bams="1.6")

load("profile.error.RData")
max.train.size <- 30
pids <- rownames(profile.error[[1]][[1]][[1]])
source("algos.in.tables.R")

error.unseen <- data.frame()
for(algo in c("flsa.norm","cghseg.k","dnacopy.sd","glad.lambdabreak")){
  for(train.size in 1:max.train.size){
    cat(sprintf("%3d / %3d train %20s\n",
                train.size,max.train.size,algo))
    nfolds <- floor(length(pids)/train.size)
    set.seed(1)
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
        train.list <- profile.error[[train.name]][[algo]]
        train.err <- colSums(with(train.list,{
          (fp+fn)[train.pids,,drop=FALSE]
        }),na.rm=TRUE)
        picked <- pick.best.index(train.err)
        for(test.name in names(profile.error)){
          test.list <- profile.error[[test.name]][[algo]]
          evec <- c()
          for(type in c("fp","fn")){
            evec[[type]] <- sum(test.list[[type]][test.pids,picked,drop=FALSE],
                                na.rm=TRUE)
          }
          evec <- c(evec,error=sum(evec))
          test.counts <-
            colSums(profile.error[[test.name]][[algo]]$anns[test.pids,])
          denoms <- c(fp=sum(test.counts[c("0breakpoints","1breakpoint")]),
                      fn=sum(test.counts[c(">0breakpoints","1breakpoint")]),
                      error=sum(test.counts))
          error.unseen <- rbind(error.unseen,{
            data.frame(algo,train.size,train.name,test.name,
                       error=evec/denoms,type=names(evec),
                       row.names=NULL)
          })
        }
      }
    }
  }
}

save(error.unseen, file="error.unseen.RData")
