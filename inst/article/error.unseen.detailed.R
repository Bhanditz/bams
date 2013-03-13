works_with_R("2.15.2",xtable="1.7.0",bams="1.6")

load("profile.error.RData")

detailed <- profile.error$detailed
pids <- rownames(detailed[[1]][[1]])

n.train.profiles <- 10
n.folds <- floor(length(pids)/n.train.profiles)

error.unseen.detailed <- data.frame()
source("algos.in.tables.R")
for(algo in algos.in.tables){
  set.seed(2)
  err.list <- detailed[[algo]]
  err.vec <- function(selection){
    with(err.list,colSums((fp+fn)[selection,,drop=FALSE]),na.rm=TRUE)
  }
  fold <- sample(rep(1:n.folds,l=length(pids)))
  ## We calculate the same thing in annotators.RData but here it is
  ## faster since we don't do it for all combinations of orig,
  ## detailed, train, test, all possible train sizes.
  for(train.fold in 1:n.folds){
    is.train <- fold == train.fold
    is.test <- !is.train
    train.err <- err.vec(is.train)
    picked <- pick.best.index(train.err)
    test.err <- err.vec(is.test)[picked]
    test.anns <- sum(err.list$anns[is.test,])
    error.unseen.detailed <- rbind(error.unseen.detailed,{
      data.frame(algo, train.fold, error=test.err/test.anns)
    })
  }
}

save(error.unseen.detailed, file="error.unseen.detailed.RData")
