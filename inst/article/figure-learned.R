works_with_R("2.15.2",neuroblastoma="1.0",bams="1.6",plyr="1.8",ggplot2="0.9.3")

source("signal.colors.R")
source("breakpoint.colors.R")
load("annotation.sets.RData")
load("zzz.stats.RData")

data(neuroblastoma)
original <- annotation.sets$original
original$annotation <- factor(original$annotation,
   c(">0breakpoints","1breakpoint","0breakpoints"))
                              

interesting <- "504"#"360"#"40"#"504"
counts <- ddply(original,.(profile.id),summarize,
      positive=sum(annotation==">0breakpoints"),
      negative=sum(annotation=="0breakpoints"))
profile.list <- with(neuroblastoma,split(profiles,profiles$profile.id))
space <- do.call(rbind,lapply(profile.list,function(df){
  chr.list <- split(df,df$chr)
  med <- median(do.call(c,lapply(chr.list,function(chr){
    diff(chr$pos)
  })))
  data.frame(profile.id=df$pro[1],
             median.spacing=med,
             probes=nrow(df))
}))
space.range <- range(space$med)
print(space.range/1e3)
probes <- sapply(profile.list,nrow)
ilist <- subset(counts,positive==3 & negative==3)$profile.id
algos <- c(cghseg.k="cghseg.k, pelt.n\nglobal model",
           ##pelt.n="pelt.n",
           dnacopy.sd="dnacopy.sd\nglobal model",
           cghseg.mBIC="cghseg.mBIC\ndefault model",
           pelt.default="pelt.default\ndefault model",
           dnacopy.default="dnacopy.default\ndefault model")
blank <- names(algos)==""
names(algos)[blank] <- algos[blank]
err.mat <- matrix(NA,length(ilist),length(algos),
                  dimnames=list(ilist,names(algos)))
for(algo in names(algos)){
  m <- all.stats[[algo]]$errors
  e <- apply(m,1,sum,na.rm=TRUE)
  min.param <- names(which.min(e))
  by.pro <- m[min.param,as.character(ilist),,drop=FALSE]
  err.mat[,algo] <- apply(by.pro,2,sum,na.rm=TRUE)
}
## Some filters on which profiles would be illustrative to
## display. i.e. these profiles are representative of the results of
## the learning algorithm.
low.density <- probes[as.character(ilist)] < 10000
cghseg.wins <- err.mat[,"cghseg.k"]<err.mat[,"dnacopy.sd"]
train.works <- err.mat[,"dnacopy.sd"]<err.mat[,"dnacopy.default"]
ilist2 <- rownames(err.mat)[cghseg.wins & low.density & train.works]
par(ask=FALSE)
plotone <- function(interesting){
  profile <- subset(profile.list[[interesting]],chromosome!="Y")
  ##  profile <- profile.list[[interesting]]
  profile$chromosome <- factor(profile$chromosome,c(1:22,"X"))
  ann <- subset(original,profile.id==interesting)
  unique.positions <-
    ddply(profile,.(chromosome,position),summarize,logratio=mean(logratio))

  smooth.list <- list()
  for(algo in names(algos)){
    m <- all.stats[[algo]]$errors
    is.train <- dimnames(m)[[2]]!=interesting
    e <- apply(m[,is.train,,drop=FALSE],1,sum,na.rm=TRUE)
    min.param <- as.numeric(names(which.min(e)))
    fun <- smoothers[[algo]]
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

  ggplot(smooth.df)+
    geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                  data=ann,alpha=1)+
    geom_vline(aes(xintercept=position/1e6),data=bkpts)+
    geom_point(aes(position/1e6,logratio),pch=1,colour="black")+
    scale_fill_manual("annotation",values=breakpoint.colors)+
    scale_x_continuous("position on chromosome (mega base pairs)",
                       breaks=c(100,200))+
    scale_y_continuous(breaks=c(-1,0,1),limits=c(-1,1))+
    facet_grid(algorithm~chromosome,scales="free_x",space="free",labeller=llab)+
    geom_line(aes(position/1e6,smooth),colour=signal.colors[["estimate"]],
              size=1.5)+
    theme_bw()+ # bioinformatics
    theme(panel.margin=unit(0,"lines"))
}

p <- plotone(interesting)
##pdf("figure-4-learned.pdf",height=6,width=13)
png("figure-4-learned.png",4200,2600,res=400,type="cairo")
print(p)
dev.off()
