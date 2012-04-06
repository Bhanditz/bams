load("zzz.stats.RData")
cids <- c("375","362")
stat.names <- c("errors","false.positive","false.negative")
ranges <- list(#flsa=10^c(-1,2),
               flsa.norm=10^c(-0.5,2),
               cghseg.k=10^c(-5,0),
               dnacopy.sd=10^c(-0.5,1.5),
               ##glad.haarseg=c(1e-60,1e-1),
               glad.lambdabreak=c(1e-1,1e5))
params <- lapply(all.stats[names(ranges)],"[[","parameters")
errs <- do.call(rbind,lapply(names(params),function(a){
  df <- do.call(rbind,lapply(stat.names,function(s){
    ar <- all.stats[[a]][[s]]
    m <- apply(ar,1:2,mean,na.rm=TRUE)
    ind <- do.call(rbind,lapply(cids,function(cid){
      data.frame(value=m[,cid],profile.id=paste("local model for profile",cid))
    }))
    glob <- data.frame(value=apply(ar,1,mean,na.rm=TRUE),
                       profile.id="global model")
    data.frame(parameter=as.numeric(params[[a]]),algorithm=a,statistic=s,
               rbind(glob,ind))
  }))
  curves <- subset(df,ranges[[a]][1] <= parameter & parameter <= ranges[[a]][2])
  transform(curves,percent=value*100)
}))
library(plyr)
library(bams)
picked <-
  ddply(subset(errs,statistic=="errors"),.(profile.id,algorithm),function(d){
  d[pick.best.index(d$value),]
})
plotted.min <- daply(errs,.(algorithm),with,min(parameter))
glob.err <- subset(errs,profile.id=="global model" & statistic=="errors")
min.labs <- subset(picked,profile.id=="global model")
min.labs$param.min <- plotted.min[as.character(min.labs$algorithm)]
min.labs$percent.text <- sprintf("%.1f",min.labs$percent)
alg.order <- as.character(with(min.labs,algorithm[order(percent)]))
errs$algorithm <- factor(errs$algorithm,alg.order)
picked$algorithm <- factor(picked$algorithm,alg.order)
min.labs$algorithm <- factor(min.labs$algorithm,alg.order)
library(ggplot2)
library(grid)
library(RColorBrewer)
cols <- brewer.pal(7,"Set1")
p <- ggplot(errs,aes(log10(parameter),percent))+
  geom_vline(aes(xintercept=log10(parameter)),colour="grey",lwd=2,
    data=subset(picked,profile.id=="global model",select=-profile.id))+
  geom_line(aes(colour=statistic,linetype=statistic),lwd=1.5)+
  geom_point(data=picked,size=4)+
  geom_segment(aes(x=log10(param.min),xend=log10(parameter),
                   y=percent,yend=percent),data=min.labs)+
  geom_text(aes(log10(param.min),percent,label=percent.text),
            data=min.labs,hjust=0,vjust=-1/2)+
  facet_grid(profile.id~algorithm,scales="free_x")+
  scale_colour_manual(values=c(#false.positive=cols[6],
                        false.positive="grey50",
                        false.negative=cols[1],
                        #errors=cols[7]
                        errors="black"
                        ))+
  scale_linetype_manual(values=c(errors=1,false.positive=2,false.negative=3))+
  xlab("log10(smoothing parameter)")+
  theme_bw()+ # bioinformatics
  opts(panel.margin=unit(0,"lines"))+
  ylab("percent incorrectly predicted annotations in training set")
pdf("figure-learning-curves.pdf",height=7.4,width=12)
print(p)
##grid.text("log10(smoothing parameter)",0.4,0,vjust=-1)
grid.text("<- more breakpoints",0.2,0,hjust=0,vjust=-1)
grid.text("fewer breakpoints ->",0.7,0,hjust=1,vjust=-1)
dev.off()
