works_with_R("2.15.2",neuroblastoma="1.0",ggplot2="0.9.3",
             bams="1.6",flsa="1.3",plyr="1.8",RColorBrewer="1.0.5")

load("annotation.sets.RData")
source("signal.colors.R")
source("breakpoint.colors.R")

data(neuroblastoma)
interesting <- "375"
profile <- subset(neuroblastoma$profiles,
                  profile.id==interesting & chromosome!="Y")
ann <- subset(annotation.sets$original, profile.id==interesting)
ann$annotation <- factor(ann$annotation,c(">0breakpoints","0breakpoints"))
unique.positions <-
  ddply(profile,.(chromosome,position),summarize,logratio=mean(logratio))
lvals <- c(0.5,7.5,10)
smooth <- ddply(unique.positions,.(chromosome),function(d){
  model <- flsa(d$logratio)
  sol <- flsaGetSolution(model,lambda2=lvals)
  do.call(rbind,lapply(seq_along(lvals),function(i){
    data.frame(d,smooth=sol[i,],lambda=lvals[i])
  }))
}) 
bkpts <- ddply(smooth,.(chromosome,lambda),function(d){
  subset(with(d,{
    data.frame(position=position[-length(position)]+diff(position)/2,
               breakpoint=diff(d$smooth))
  }),breakpoint!=0)
})

library(grid)
llab <- function(var,val){
  if(var=="lambda")sprintf("lambda = %s",val)
  else as.character(val)
}

p <- ggplot(smooth)+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                data=ann,alpha=1)+
  geom_vline(aes(xintercept=position/1e6),colour="black",data=bkpts)+
  geom_point(aes(position/1e6,logratio),pch=1,colour="black")+
  scale_fill_manual("annotation",values=breakpoint.colors)+
  scale_x_continuous("position on chromosome (mega base pairs)",
                     breaks=c(100,200))+
  scale_y_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
  facet_grid(lambda~chromosome,scales="free_x",space="free",labeller=llab)+
  geom_line(aes(position/1e6,smooth),colour=signal.colors[["estimate"]],
            size=1.5)+
  theme_bw()+ # bioinformatics
  theme(panel.margin=unit(0,"lines"))

pdf("figure-1-smoothing.pdf",height=5.3,width=13)
print(p)
dev.off()
