library(bams)
load("zzz.stats.RData")
results <- data.frame()
display.names <- c(cghseg.k="cghseg.k, pelt.n",
                   "flsa.norm","dnacopy.sd","glad.lambdabreak")
blank <- names(display.names)==""
names(display.names)[blank] <- display.names[blank]
for(algorithm in names(display.names)){
  stat <- all.stats[[algorithm]]
  for(N in c(1:30)){
    set.seed(2)
    est <- (1-unseen.profile.error(stat,N)) * 100
    results <- rbind(results,data.frame(statistic=colnames(est),
                                        mean=apply(est,2,mean),
                                        sd=apply(est,2,sd),
                                        training.set.profiles=N,
                                        algorithm=display.names[[algorithm]]))
  }
}
library(ggplot2)
source("algo.colors.R")
library(grid)
algo.labels <- data.frame(training.set.profiles=29.5,
                          mean=c(82.5,89.5,92.5,99),
  algorithm=c("glad.lambdabreak","dnacopy.sd",
    "flsa.norm","cghseg.k, pelt.n"))
p <- ggplot(subset(results,statistic=="errors"),
            aes(training.set.profiles,mean))+
  geom_vline(x=10)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,group=algorithm,fill=algorithm),
              alpha=1/4)+
  geom_line(aes(group=algorithm,colour=algorithm),lwd=1.5)+
  scale_x_continuous("Annotated profiles in global model training set",
                     breaks=c(1,5,10,15,20,25,30))+
  scale_y_continuous(paste("Percent of correctly predicted annotations",
                           "on test set profiles"),
                     breaks=seq(80,100,by=2))+
  scale_fill_manual(values=algo.colors,guide="none")+
  scale_colour_manual(values=algo.colors,guide="none")+
  coord_cartesian(xlim=c(1,30),ylim=c(80,100))+
  theme_bw()+
  geom_text(aes(colour=algorithm,label=algorithm),data=algo.labels,
            hjust=1)
pdf("figure-5-kinetics.pdf",width=5.1,height=5.2)
print(p)
dev.off()
