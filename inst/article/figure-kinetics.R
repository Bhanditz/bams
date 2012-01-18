library(bams)
load("zzz.stats.RData")
results <- data.frame()
for(algorithm in c("cghseg.k","flsa.norm","dnacopy.sd","glad.lambdabreak")){
  stat <- all.stats[[algorithm]]
  for(N in c(1:30)){
    set.seed(2)
    est <- unseen.profile.error(stat,N) * 100
    results <- rbind(results,data.frame(statistic=colnames(est),
                                        mean=apply(est,2,mean),
                                        sd=apply(est,2,sd),
                                        training.set.profiles=N,
                                        algorithm))
  }
}
library(ggplot2)
algo.labels <- data.frame(training.set.profiles=29.5,
                          mean=c(17,10,7,1),
  algorithm=c("glad.lambdabreak","dnacopy.sd",
    "flsa.norm","cghseg.k"))
p <- ggplot(subset(results,statistic=="errors"),
            aes(training.set.profiles,mean))+
  geom_vline(x=10)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,group=algorithm,fill=algorithm),
              alpha=1/4)+
  geom_line(aes(group=algorithm,colour=algorithm),lwd=1.5)+
  scale_x_continuous("Annotated profiles in global model training set",
                     breaks=c(1,5,10,15,20,25,30))+
  scale_y_continuous(paste("Percent of incorrectly predicted annotations",
                           "on test set profiles"),
                     breaks=seq(0,20,by=2))+
  scale_fill_discrete(legend=FALSE)+
  scale_colour_discrete(legend=FALSE)+
  coord_cartesian(xlim=c(1,30),ylim=c(0,20))+
  theme_bw()+
  opts(axis.title.x=theme_text(vjust = 0))+
  geom_text(aes(colour=algorithm,label=algorithm),data=algo.labels,
            hjust=1)
pdf("figure-kinetics.pdf")
print(p)
dev.off()
