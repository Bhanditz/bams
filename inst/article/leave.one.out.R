works_with_R("2.15.2",bams="1.6")

load("annotation.sets.RData")
load("fp.fn.list.RData")

## Calculate some error estimates when we just leave one annotation in
## each profile out, then use local or global models.

pids <- levels(annotation.sets[[1]]$profile.id)
set.seed(1)
nboot <- 10
leave.one.out <- data.frame()
for(set.name in names(annotation.sets)){
  set <- annotation.sets[[set.name]]
  ann.list <- split(set,set$profile.id)
  algo.list <- fp.fn.list[[set.name]]
  for(algo in names(algo.list)){
    cat(sprintf("%10s %20s\n",set.name,algo))
    err.list <- list()
    for(type in c("fp","fn")){
      err.list[[type]] <-
        lapply(split(data.frame(algo.list[[algo]][[type]]),set$profile.id),
               as.matrix)
    }
    for(i in 1:nboot){
      ## errors pids x parameters, with held out data.
      nparam <- ncol(err.list[[1]][[1]])
      train.error.mat <- matrix(NA,length(pids),nparam)
      rownames(train.error.mat) <- pids
      test.list <- list(anns=rep(NA,length(pids)))
      names(test.list$anns) <- pids
      for(type in c("fp","fn")){
        test.list[[type]] <- train.error.mat #shortcut.
      }
      ## First pass through profiles: hold out one annotation and
      ## record what is the train and test error for each model.
      for(pid in pids){
        anns <- ann.list[[pid]]
        err <- lapply(err.list,"[[",pid)
        test <- sample(nrow(anns),1)
        test.list$anns[pid] <- if(any(is.na(err[[1]]))){
          NA
        }else{
          anns[test,"annotation"]
        }
        for(type in c("fp","fn")){
          test.list[[type]][pid,] <- err[[type]][test,]
        }
        train.error.mat[pid,] <- with(err,colSums((fp+fn)[-test,,drop=FALSE]))
      }
      ## Train a global model.
      picked <- pick.best.index(colSums(train.error.mat,na.rm=TRUE))
      test.list$global <- sapply(c("fp","fn"),function(type){
        sum(test.list[[type]][,picked],na.rm=TRUE)
      })
      ## Second pass: train a local model and record errors.
      test.list$local <- rowSums(sapply(pids,function(pid){
        loc.err <- train.error.mat[pid,]
        if(any(is.na(loc.err))){
          rep(NA,2)
        }else{
          picked <- pick.best.index(loc.err)
          sapply(c("fp","fn"),function(type){
            t(test.list[[type]][pid,picked])
          })
        }
      }),na.rm=TRUE)
      test.counts <- table(factor(test.list$anns,
                    c("0breakpoints",">0breakpoints","1breakpoint")))
      denoms <- c(fp=sum(test.counts[c("0breakpoints","1breakpoint")]),
                  fn=sum(test.counts[c(">0breakpoints","1breakpoint")]),
                  error=sum(test.counts))
      for(model in c("global","local")){
        error <- test.list[[model]]
        error <- c(error,error=sum(error))
        leave.one.out <- rbind(leave.one.out,{
          data.frame(set.name,algo,i,model,
                     error=error/denoms,
                     type=names(error))
        })
      }
    }
  }
}

save(leave.one.out, file="leave.one.out.RData")
                 
