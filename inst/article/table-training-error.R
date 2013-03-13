works_with_R("2.15.2",bams="1.6",reshape2="1.2.2",xtable="1.7.0",plyr="1.8")

load("profile.error.RData")
source("squares.R")

local.training.error <- data.frame()
for(set.name in names(profile.error)){
  set <- profile.error[[set.name]]
  for(algo in names(set)){
    L <- set[[algo]]
    pick.vec <- with(L, apply(fp+fn, 1, pick.best.index))
    select.mat <- cbind(seq_along(pick.vec), pick.vec)
    err <- sapply(c("fp","fn"),function(type){
      sum(L[[type]][select.mat],na.rm=TRUE)
    })
    err <- c(err,error=sum(err))
    counts <- colSums(L$ann)
    denoms <- c(fp=sum(counts[c("0breakpoints","1breakpoint")]),
                fn=sum(counts[c(">0breakpoints","1breakpoint")]),
                error=sum(counts))
    local.training.error <- rbind(local.training.error,{
      data.frame(set.name, algo, type=names(err), error=err/denoms*100,
                 row.names=NULL)
    })
  }
}

## save the percents for textual discussion.
sets <- ddply(subset(local.training.error,type=="error"),.(set.name),
              function(df)df[which.min(df$err),])
for(i in 1:nrow(sets)){
  r <- sets[i,]
  f <- sprintf("local-%s.tex", r$set)
  cat(round(r$error, 1), file=f)
}

wide <- dcast(local.training.error, algo~set.name+type)
wide$algo <- squares(as.character(wide$algo))

xt <- xtable(wide[order(wide[,"detailed_error"]),],digits=1,
             align="rrrrrrrr")
colnames(xt) <- c("algorithm",rep(c("error","FP","FN"),2))
tex <- print(xt,include.rownames=FALSE,floating=FALSE,
             sanitize.text.function=identity)
tex <- sub("hline","hline & \\multicolumn{3}{c}{original} & \\multicolumn{3}{c}{detailed} \\\\",tex,fixed=TRUE)
cat(tex,file="table-training-error.tex")
