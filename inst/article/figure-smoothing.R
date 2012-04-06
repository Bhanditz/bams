data(neuroblastoma,package="neuroblastoma")
interesting <- "375"
profile <- subset(neuroblastoma$profiles,profile.id==interesting)
ann <- subset(neuroblastoma$annotations,profile.id==interesting)
library(bams)
library(flsa)
library(plyr)
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
library(ggplot2)
library(grid)
llab <- function(var,val){
  if(var=="lambda")sprintf("lambda = %s",val)
  else as.character(val)
}
library(RColorBrewer)
cols <- brewer.pal(7,"Set1")
smooth.col <- cols[2]
smooth.col <- "#00aaff"
p <- ggplot(smooth)+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                data=ann,alpha=1)+
  geom_point(aes(position/1e6,logratio),pch=1,colour="black")+
  scale_fill_manual("annotation",values=c(breakpoint=cols[1],
                                   #normal=cols[6]
                                   normal="grey50"
                                   ))+
  scale_x_continuous("position on chromosome (mega base pairs)",
                     breaks=c(100,200))+
  scale_y_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
  facet_grid(lambda~chromosome,scales="free_x",space="free",labeller=llab)+
  geom_line(aes(position/1e6,smooth),colour=smooth.col,size=1.5)+
  geom_vline(aes(xintercept=position/1e6),size=1,
             colour=smooth.col,linetype="dashed",data=bkpts)+
  theme_bw()+ # bioinformatics
  opts(panel.margin=unit(0,"lines"))
pdf("figure-smoothing.pdf",height=4.1,width=12)
print(p)
dev.off()
