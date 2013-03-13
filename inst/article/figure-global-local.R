works_with_R("2.15.2",ggplot2="0.9.3",plyr="1.8")

load("leave.one.out.RData")
source("algo.colors.R")
algo.colors[algo.colors=="#FFFFFF"] <- NA
all.err <- ddply(leave.one.out,.(set.name,algo,model,type),summarize,
                 mean=mean(error),sd=sd(error))

error <- subset(all.err,type=="error")
glob.orig <- subset(error,model=="global" & set.name=="original")
error$algo <- factor(error$algo,{
  with(glob.orig,as.character(algo)[order(mean,decreasing=TRUE)])
})
error$set.name <- factor(error$set.name,c("detailed","original"))
library(grid)
squares <- data.frame(algo=levels(error$algo),mean=0,set.name="detailed")
p <- ggplot(error,aes(mean,algo,colour=model))+
  geom_point(size=4,pch="|")+
  geom_point(aes(fill=algo,colour=NULL), data=squares, colour=NA,
             shape=22, size=5)+
  geom_segment(aes(mean-sd,xend=mean+sd,yend=algo,size=model))+
  facet_grid(.~set.name,labeller=function(var,val){
    sprintf("annotations: %s",val)
  })+
  scale_size_manual(values=c(global=2,local=1))+
  scale_colour_manual(values=c(global="black",local="red"))+
  scale_fill_manual(values=algo.colors, guide="none")+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  ylab("algorithm")+
  ##coord_cartesian(xlim=c(0,0.85))+
  xlab(expression(
      "Leave-one-out test error (mean " %+-% "1 standard deviation)"))

sets <- ddply(error,.(set.name),
              function(df)df[which.min(df$mean),])
for(i in 1:nrow(sets)){
  r <- sets[i,]
  f <- sprintf("global-test-%s.tex", r$set)
  cat(round(r$mean*100, 1), file=f)
}

pdf("figure-global-local.pdf",w=8,h=4)
print(p)
dev.off()
