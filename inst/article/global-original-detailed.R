works_with_R("2.15.2",plyr="1.8")

load("global.error.RData")

sets <- ddply(global.error, .(set), function(df){
  df[which.min(df$error),]
})

for(i in 1:nrow(sets)){
  r <- sets[i,]
  f <- sprintf("global-%s.tex", r$set)
  cat(round(r$error*100, 1), file=f)
}

