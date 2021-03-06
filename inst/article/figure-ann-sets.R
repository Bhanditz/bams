works_with_R("2.15.2",bams="1.6",ggplot2="0.9.3",neuroblastoma="1.0")
load("annotation.sets.RData")
source("breakpoint.colors.R")

data(neuroblastoma)
pid <- "8"

ann.df <- data.frame()
for(set.name in names(annotation.sets)){
  set <- annotation.sets[[set.name]]
  these <- subset(set,profile.id%in%pid)
  these$set <- set.name
  ann.df <- rbind(ann.df,these)
}
ann.df$annotation <- factor(ann.df$annotation,
   c(">0breakpoints","1breakpoint","0breakpoints"))

profile <- subset(neuroblastoma$profiles,
                  profile.id==pid & chromosome%in%unique(ann.df$chr))

lfun <- function(var,val){
  if(var=="set"){
    sprintf("annotations: %s",val)
  }else as.character(val)
}
library(grid)
p <- ggplot()+
  geom_blank(aes(min/1e6),y=0,data=ann.df)+
  geom_blank(aes(max/1e6),y=0,data=ann.df)+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),data=ann.df)+
  geom_point(aes(position/1e6,logratio),data=profile)+
  facet_grid(set~chromosome,scales="free",space="free_x",
             labeller=lfun)+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  scale_x_continuous("position on chromosome (mega base pairs)",
                     breaks=c(100,200))+
  scale_fill_manual(values=breakpoint.colors)

png("figure-ann-sets.png",2800,1600,res=400,type="cairo")
print(p)
dev.off()
