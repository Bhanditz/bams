works_with_R("2.15.2",ggplot2="0.9.3",plyr="1.8")

load("error.unseen.RData")
source("algo.colors.R")

err.only <- subset(error.unseen, type == "error")
stats <- ddply(err.only,.(train.name,test.name,algo,train.size),summarize,
               mean=(1-mean(error))*100,sd=sd(error)*100)

library(grid)
l <- function(mean,algorithm,algo=algorithm,test.name="original",
              train.name="original",train.size=29.5){
  data.frame(mean,algorithm,test.name,train.name,train.size,algo)
}
algo.labels <- rbind(l(82.5,"glad.lambdabreak"),
                     l(90,"dnacopy.sd"),
                     l(92.5,"flsa.norm"),
                     l(92.5,"gada",train.size=20),
                     l(96,"cghseg.k, pelt.n","cghseg.k",
                       test.name="detailed"))
vline <- data.frame(train.size=10,test.name="detailed",train.name="detailed")
toplot <- subset(stats, algo%in%as.character(algo.labels$algo))
p <- ggplot(toplot,aes(train.size,mean))+
  geom_vline(aes(xintercept=train.size),data=vline)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,group=algo,fill=algo),
              alpha=1/4)+
  geom_line(aes(group=algo,colour=algo),lwd=1.5)+
  scale_x_continuous("Annotated profiles in global model training set",
                     breaks=c(1,5,10,15,20,25))+
  scale_y_continuous(paste("Percent of correctly predicted annotations",
                           "on test set profiles"))+
  scale_fill_manual(values=algo.colors,guide="none")+
  scale_colour_manual(values=algo.colors,guide="none")+
  coord_cartesian(xlim=c(1,30),ylim=c(78,100))+
  theme_bw()+
  geom_text(aes(colour=algorithm,label=algorithm),data=algo.labels,
            hjust=1)+
  theme(panel.margin=unit(0,"cm"))+
  facet_grid(test.name~train.name,labeller=function(var,val){
    sprintf("%s: %s",sub(".name","",var),val)
  })

pdf("figure-test-error-train-profiles.pdf",h=5)
print(p)
dev.off()
