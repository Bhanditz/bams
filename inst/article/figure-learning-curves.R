works_with_R("2.15.2",bams="1.6",plyr="1.8")

source("breakpoint.colors.R")
load("profile.error.RData")
load("parameter.list.RData")

fp.fn.colors <- c(false.positive="skyblue",
                  false.negative="#E41A1C",
                  errors="black",
                  imprecision="black")

set.name <- "original"
loc <- profile.error[[set.name]]

ranges <- list(#flsa=10^c(-1,2),
               flsa.norm=10^c(-0.3,2),
               cghseg.k=10^c(-5,-0.5),
               ##glad.lambdabreak=c(1e-1,1e5),
               dnacopy.sd=10^c(-0.5,1.5))
               ##glad.haarseg=c(1e-60,1e-1),
cids <- c("375","362")
display.name <- c(cghseg.k="cghseg.k, pelt.n")
##load("zzz.stats.RData")
##stat.names <- c("errors","false.positive","false.negative")
## params <- lapply(all.stats[names(ranges)],"[[","parameters")
## errs <- do.call(rbind,lapply(names(params),function(a){
##   df <- do.call(rbind,lapply(stat.names,function(s){
##     ar <- all.stats[[a]][[s]]
##     m <- apply(ar,1:2,mean,na.rm=TRUE)
##     ind <- do.call(rbind,lapply(cids,function(cid){
##       data.frame(value=m[,cid],profile.id=paste("local model for profile",cid))
##     }))
##     glob <- data.frame(value=apply(ar,1,mean,na.rm=TRUE),
##                        profile.id="global model")
##     algorithm <- if(a %in% names(display.name)){
##       display.name[[a]]
##     }else{
##       a
##     }
##     data.frame(parameter=as.numeric(params[[a]]),
##                algorithm,
##                statistic=s,
##                rbind(glob,ind))
##   }))
##   curves <- subset(df,ranges[[a]][1] <= parameter & parameter <= ranges[[a]][2])
##   transform(curves,percent=value*100)
## }))
errs <- data.frame()
for(a in names(ranges)){
  algorithm <- if(a %in% names(display.name)){
    display.name[[a]]
  }else{
    a
  }
  info <- loc[[a]]
  parameter <- parameter.list[[a]]
  r <- ranges[[a]]
  show.parameter <- r[1] < parameter & parameter < r[2]
  parameter <- parameter[show.parameter]
  selectors <- list("global model"=rownames(info$tp))
  for(cid in cids){
    mname <- sprintf("local model for profile %s",cid)
    selectors[[mname]] <- cid
  }
  stat.funs <- list(errors=function(L)L$fp+L$fn,
                    false.positive=function(L)L$fp,
                    false.negative=function(L)L$fn)
  for(statistic in names(stat.funs)){
    fun <- stat.funs[[statistic]]
    err.mat <- fun(info)[,show.parameter]
    for(profile.id in names(selectors)){
      selected <- selectors[[profile.id]]
      sub.mat <- err.mat[selected,,drop=FALSE]
      n.anns <- sum(info$anns[selected,])
      percent <- colSums(sub.mat)/n.anns * 100
      errs <- rbind(errs,{
        data.frame(parameter,algorithm,statistic,percent,profile.id)
      })
    }
  }
}
picked <-
  ddply(subset(errs,statistic=="errors"),.(profile.id,algorithm),function(d){
  d[pick.best.index(d$percent),]
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
  geom_line(aes(colour=statistic,size=statistic))+
  geom_point(data=picked,size=4)+
  geom_segment(aes(x=log10(param.min),xend=log10(parameter),
                   y=percent,yend=percent),data=min.labs)+
  geom_text(aes(log10(param.min),percent,label=percent.text),
            data=min.labs,hjust=0,vjust=-1/2)+
  facet_grid(profile.id~algorithm,scales="free_x")+
  scale_colour_manual(values=fp.fn.colors)+
  scale_size_manual(values=c(
                      false.positive=1.5,
                      false.negative=0.8,
                      errors=2.5))+
  ##scale_linetype_manual(values=c(errors=1,false.positive=2,false.negative=3))+
  xlab("log10(smoothing parameter)")+
  theme_bw()+ # bioinformatics
  theme(panel.margin=unit(0,"lines"))+
  ylab("percent incorrectly predicted annotations in training set")
  

pdf("figure-2-learning-curves.pdf",height=6,width=10)
print(p)
##grid.text("log10(smoothing parameter)",0.4,0,vjust=-1)
grid.text("<- more breakpoints",0.1,0,hjust=0,vjust=-1)
grid.text("fewer breakpoints ->",0.825,0,hjust=1,vjust=-1)
dev.off()
