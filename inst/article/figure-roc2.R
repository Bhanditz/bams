works_with_R("2.15.2",ggplot2="0.9.3",plyr="1.8",bams="1.6")

load("global.error.RData")
source("algo.class.R")
source("algo.colors.R")

global.error$class <- algo.class[as.character(global.error$algo)]
na.class <- subset(global.error,is.na(class))
if(nrow(na.class)){
  print(na.class)
  stop("some NA class values")
}

dots <- ddply(global.error,.(set,algo),function(df){
  df[pick.best.index(df$error),]
})
dots$model <- ifelse(is.na(dots$parameter),"default","global")

l <- function(algo,set,fpr,tpr,label=algo){
  data.frame(algo,set,class=algo.class[algo],fpr,tpr,label)
}
dlabs <- rbind(l("cghseg.k","detailed",0.07,0.98),
               l("pelt.n","detailed",0.1,0.85),
               l("flsa.norm","detailed",0.22,0.65),
               l("flsa","original",0.12,0.67),
               l("pelt.default","original",0.15,0.58),
               l("cghseg.mBIC","detailed",0.4,0.87),
               l("dnacopy.alpha","detailed",0.18,0.97),
               l("dnacopy.default","detailed",0.4,0.88,"dnacopy\ndefault"),
               l("dnacopy.prune","detailed",0.37,0.78),
               l("dnacopy.sd","original",0.18,0.68),
               l("gada","original",0.07,0.83),
               l("glad.haarseg","detailed",0.1,0.95),
               l("glad.lambdabreak","original",0.26,0.75),
               l("glad.default","detailed",0.38,0.92),
               l("glad.MinBkpWeight","original",0.3,0.87))
algo.colors <- sub("FFFFFF","000000",algo.colors)
library(grid)
p <- ggplot(global.error,aes(fpr,tpr))+
  facet_grid(set~class,labeller=function(var,val){
    if(var=="set"){
      sprintf("annotations: %s",val)
    }else{
      as.character(val)
    }
  })+
  geom_path(aes(colour=algo),lwd=1.5)+
  geom_path(aes(group=algo),lty="dashed")+
  geom_point(aes(shape=model),data=dots,size=4)+
  geom_text(aes(colour=algo,label=label),data=dlabs,size=3)+
  coord_cartesian(xlim=c(0,0.5),ylim=c(0.5,1))+
  scale_x_continuous(paste("False positive rate = ",
            "probability(predict too many breakpoints | 0 or 1 breakpoint)",
                           sep=""),
                     breaks=seq(0,0.4,by=0.1))+
  scale_y_continuous(paste("True positive rate =\n",
           "probability(predict breakpoint | at least 1 breakpoint)",
                           sep=""),
                     breaks=seq(0.6,1,by=0.1))+
  theme_bw()+
  scale_colour_manual(values=algo.colors,guide="none")+
  theme(panel.margin=unit(0,"cm"))+
  scale_shape_manual(values=c(default=19,global=1))

pdf("figure-roc2.pdf",w=8.5,h=5)
print(p)
dev.off()
