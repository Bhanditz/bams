works_with_R("2.15.2",ggplot2="0.9.3",neuroblastoma="1.0",bams="1.6")

load("profile.errors.RData")
load("segmentation.list.RData")
load("annotation.sets.RData")
source("breakpoint.colors.R")
source("signal.colors.R")

## show weird profiles for which flsa does better than cghseg.k.
e.list <- profile.error$detailed
fn <- e.list$cghseg.k$fn
pids <- rownames(fn)[apply(fn,1,min)>0]

data(neuroblastoma)
p.list <- split(neuroblastoma$pro, neuroblastoma$pro$pro)
a.list <- split(annotation.sets$detailed, annotation.sets$detailed$pro)

p.df <- do.call(rbind, p.list[pids])
a.df <- do.call(rbind, a.list[pids])
brk.df <- data.frame()
for(pid in pids){
  for(algo in c("cghseg.k","flsa")){
    e <- with(e.list[[algo]],(fp+fn)[pid,])
    picked <- pick.best.index(e)
    chr.list <- segmentation.list[[pid]][[algo]]$breakpoints
    for(chromosome in names(chr.list)){
      bvec <- chr.list[[chromosome]][[picked]]
      if(length(bvec)){
        brk.df <- rbind(brk.df,{
          data.frame(profile.id=pid, algo, chromosome, position=bvec)
        })
      }
    }
  }
}
library(grid)

for(algorithm in unique(brk.df$algo)){
p <- ggplot()+
  geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6, fill=annotation), data=a.df)+
  scale_fill_manual(values=breakpoint.colors)+
  geom_point(aes(position/1e6, logratio), data=p.df)+
  geom_vline(aes(xintercept=position/1e6),
             data=subset(brk.df,algo==algorithm),
             colour=signal.colors[["estimate"]], linetype="dashed")+
  facet_grid(profile.id~chromosome,scales="free_x",space="free")+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))
png(sprintf("figure-local-breaks-%s.png",algorithm),
    4200,2600,res=400,type="cairo")
print(p)
dev.off()
}
