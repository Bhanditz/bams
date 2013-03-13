load("profile.error.RData")
load("parameter.list.RData")

## list[[set]][[algo]] is a list: fp and fn are nprofile x nparameter
## matrices, and anns are annotation counts.

global.error <- data.frame()
for(set.name in names(profile.error)){
  set <- profile.error[[set.name]]
  for(algo in names(set)){
    stats <- set[[algo]]
    parameter <- parameter.list[[algo]]
    incorrect <- with(stats,{
      colSums(fp)+colSums(fn)
    })
    error <- incorrect/sum(stats$anns)
    tpr <- with(stats,{
      colSums(tp)/sum(anns[,c("1breakpoint",">0breakpoints")])
    })
    fpr <- with(stats,{
      colSums(fp)/sum(anns[,c("1breakpoint","0breakpoints")])
    })
    global.error <- rbind(global.error,{
      data.frame(set=set.name,algo,parameter,error,tpr,fpr)
    })
  }
}

save(global.error, file="global.error.RData")
