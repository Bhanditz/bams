works_with_R("2.15.2",neuroblastoma="1.0",bams="1.6",plyr="1.8",ggplot2="0.9.3")

source("signal.colors.R")
source("breakpoint.colors.R")
load("annotation.sets.RData")
load("profile.error.RData")

data(neuroblastoma)
set.name <- "detailed"
algo.err.list <- profile.error[[set.name]]
detailed <- annotation.sets[[set.name]]
detailed$annotation <- factor(detailed$annotation,
   c(">0breakpoints","1breakpoint","0breakpoints"))

counts <- daply(detailed,.(profile.id),with,table(annotation))
profile.list <- with(neuroblastoma,split(profiles,profiles$profile.id))

probes <- sapply(profile.list,nrow)
ilist <- c(rownames(counts)[apply(counts>2,1,all)],
           rownames(counts)[probes>1e4 & counts[,"1breakpoint"]>0],
           "504")
algos <- c(cghseg.k="cghseg.k, pelt.n\nglobal model",
           dnacopy.sd="dnacopy.sd\nglobal model",
           cghseg.mBIC="cghseg.mBIC\ndefault model",
           pelt.default="pelt.default\ndefault model",
           dnacopy.default="dnacopy.default\ndefault model")

plotone <- function(interesting){
  profile <- subset(profile.list[[interesting]],chromosome!="Y")
  ##  profile <- profile.list[[interesting]]
  profile$chromosome <- factor(profile$chromosome,c(1:22,"X"))
  ann <- subset(detailed,profile.id==interesting)
  unique.positions <-
    ddply(profile,.(chromosome,position),summarize,logratio=mean(logratio))

  smooth.list <- list()
  for(algo in names(algos)){
    m <- with(algo.err.list[[algo]],fp+fn)
    is.train <- !rownames(m)%in%ilist
    e <- colSums(m[is.train,,drop=FALSE],na.rm=TRUE)
    picked <- pick.best.index(e)
    fun <- smoothers[[algo]]
    min.param <- eval(formals(fun)[[2]])[picked]
    unique.positions$smooth <- t(fun(unique.positions,min.param))
    unique.positions$algorithm <- algos[[algo]]
    smooth.list[[algo]] <- unique.positions
  }
  smooth.df <- do.call(rbind,smooth.list)
  smooth.df$algorithm <- factor(smooth.df$algorithm,algos)
  ## smooth has columns chromosome position logratio smooth lambda
  bkpts <- ddply(smooth.df,.(chromosome,algorithm),function(d){
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

  makegg <- function(chr=NULL){
    if(!is.null(chr)){
      smooth.df <- subset(smooth.df, chromosome%in%chr)
      bkpts <- subset(bkpts, chromosome%in%chr)
    }
    ggplot(smooth.df)+
       geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                     data=ann,alpha=1)+
       geom_vline(aes(xintercept=position/1e6),data=bkpts)+
       geom_point(aes(position/1e6,logratio),pch=1,colour="black")+
       scale_fill_manual("annotation",values=breakpoint.colors)+
       scale_x_continuous("position on chromosome (mega base pairs)",
                          breaks=c(100,200))+
       facet_grid(algorithm~chromosome,
                  scales="free_x",space="free",labeller=llab)+
       geom_line(aes(position/1e6,smooth),colour=signal.colors[["estimate"]],
                 size=1.5)+
       theme_bw()+ # bioinformatics
       theme(panel.margin=unit(0,"lines"))
  }
  list(all=makegg(),anns=makegg(unique(ann$chr)))
}

par(ask=FALSE)
for(pid in ilist){
  gglist <- plotone(pid)
  for(N in names(gglist)){
    p <- gglist[[N]]
    f <- sprintf("figure-supplementary-%s-%s.png",pid,N)
    cat(f,"\n")
    png(f,min(2*probes[pid],32000),2600,res=400,type="cairo")
    print(p)
    dev.off()
  }
}
