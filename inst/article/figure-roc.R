load("zzz.stats.RData")
library(bams)
roc <- data.frame()
for(algorithm in names(all.stats)){
  stat <- all.stats[[algorithm]]
  normal.anns <- sum(stat$normal.anns)
  breakpoint.anns <- sum(stat$breakpoint.anns)
  correct <- stat$errors==0
  TPR <- apply(correct,1,function(correct.mat){
    sum(correct.mat & stat$breakpoint.anns) / breakpoint.anns
  })
  FPR <- apply(stat$errors,1,function(error.mat){
    sum(error.mat & stat$normal.anns) / normal.anns
  })
  errors <- apply(stat$errors,1,sum,na.rm=TRUE)
  class <- gsub("[.].*","",algorithm)
  class <- switch(class,
                  cghFLasso="optimization",
                  flsa="optimization",
                  cghseg="optimization",
                  class)
  class <- gsub("optimization","optimization-based models",class)
  best.i <- pick.best.index(errors)
  best <- rep(FALSE,length(errors))
  best[best.i] <- TRUE
  newroc <- data.frame(class,algorithm,parameter=stat$parameters,
                       TPR,FPR,errors,best)
  roc <- rbind(roc,newroc)
}
library(lattice)
##xyplot(TPR~FPR|algorithm,roc,type="o")
no.tuning <- names(all.stats)[sapply(all.stats,function(L)dim(L$errors)[1])==1]
curves <- subset(roc,!algorithm%in%no.tuning)
## Need to reorder factor to get different looking colors in the same
## panel
algo.class <- unique(curves[,c("algorithm","class")])
algo.class <- algo.class[order(algo.class$class),]
lev.list <- lapply(levels(algo.class$class),function(lev){
  as.character(subset(algo.class,class==lev)$algorithm)
})
levs <- c()
for(i in 1:max(sapply(lev.list,length))){
  for(v in lev.list){
    if(length(v)>=i)levs <- c(levs,v[i])
  }
}
curves$algorithm <- factor(curves$algorithm,levs)
dots <- subset(roc,algorithm%in%no.tuning)
dots$vjust <- 2
dots$hjust <- -0.1
dots$label <- as.character(dots$algorithm)
colordots <- subset(curves,best)
colordots$hjust <- -0.1
colordots$vjust <- 1.1
colordots$label <- as.character(colordots$algorithm)
change <- list(hjust=c(dnacopy.alpha=1.1,glad.haarseg=0.5,flsa=-0.5,
                 glad.MinBkpWeight=0.5,dnacopy.default=0.9,
                 cghseg.mBIC=1,cghseg.k=0.05),
               vjust=c(dnacopy.alpha=0,glad.haarseg=-3,
                 glad.MinBkpWeight=5,dnacopy.default=1.7,cghseg.k=-4.7),
               label=c(flsa.norm=" flsa\nnorm",
                 dnacopy.default=" dnacopy\ndefault"))
for(coln in names(change)){
  v <- change[[coln]]
  for(N in names(v)){
    colordots[colordots$algorithm==N,coln] <- v[N]
    dots[dots$algorithm==N,coln] <- v[N]
  }
}

source("algo.colors.R")
library(ggplot2)
library(grid)
dotsize <- 6
text.cex <- 4
p <- ggplot(roc,aes(FPR,TPR))+
  facet_grid(.~class)+
  geom_path(aes(colour=algorithm),data=curves,lwd=1.5)+
  geom_path(aes(group=algorithm),data=curves,lty="dashed")+
  geom_point(data=dots,fill="black",colour="black",pch=21,size=dotsize)+
  geom_point(fill=NA,pch=21,data=colordots,size=dotsize)+
  geom_text(aes(colour=algorithm,label=label,hjust=hjust,vjust=vjust),
            data=colordots,cex=text.cex)+
  geom_text(aes(label=label,hjust=hjust,vjust=vjust),data=dots,
            cex=text.cex)+
  coord_cartesian(xlim=c(0,0.5),ylim=c(0.5,1))+
  scale_x_continuous("False positive rate = probability(predict breakpoint | normal)",
                     breaks=seq(0,0.4,by=0.1))+
  scale_y_continuous(paste("True positive rate =\n",
                           "probability(predict breakpoint | breakpoint)",
                           sep=""),
                     breaks=seq(0.5,1,by=0.1))+
  theme_bw()+
  scale_colour_manual(values=algo.colors,guide="none")
pdf("figure-roc.pdf",width=10,height=4)
print(p)
dev.off()
