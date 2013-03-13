works_with_R("2.15.2",SegAnnot="1.1")

data(profiles)

pro <- profiles$hi
chr <- pro$pro
##results <- with(chr,run.cghseg(logratio,position,13))
anns <- pro$ann
win <- function(min,max)data.frame(min=min*1e5,max=max*1e5)
windows <- rbind(##win(  65,  71),
                 win( 148, 171),
                 win( 355, 361),
                 win(1060,1065))
anns$mid <- with(anns,(min+max)/2)
chr.df <- data.frame()
ann.df <- data.frame()
##seg.df <- data.frame()
##breaks <- data.frame()
##showSegs <- c(1,2,3,5,7,13)
for(i in 1:nrow(windows)){
  w <- windows[i,]
  sm.df <- subset(results$seg,!(first.base>w$max | last.base<w$min))
  sm.df$first.base[sm.df$first.base < w$min] <- w$min
  sm.df$last.base[sm.df$last.base > w$max] <- w$max
  ##breaks <- rbind(breaks,{
    ##data.frame(subset(results$break.df,w$min < base & base < w$max &
      ##                segments%in%showSegs),i)
  ##  })
  ##seg.df <- rbind(seg.df,{
    ##data.frame(subset(sm.df,segments%in%showSegs),i)
##  })
  chr.df <- rbind(chr.df,{
    data.frame(subset(chr,w$min < position & position < w$max),i)
  })
  ann.df <- rbind(ann.df,{
    data.frame(subset(anns,w$min < mid & mid < w$max),i)
  })
}
labeller <- function(var,val){
  if(var=="segments"){
    s <- ifelse(val==1,"","s")
    sprintf("%s segment%s",val,s)
  }else{
    sprintf("window %d",val)
  }
}
br <- c(6.5,7.0,seq(15,17,by=0.5),35.5,36,106,106.5)
names(br) <- as.character(br)
p <- ggplot()+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                colour="black",data=ann.df)+
  geom_point(aes(position/1e6,logratio),data=chr.df)+
  ## geom_segment(aes(first.base/1e6,mean,xend=last.base/1e6,yend=mean),
  ##              data=seg.df,colour=signal.colors[["estimate"]],lwd=1)+
  ## geom_vline(aes(xintercept=base/1e6),data=breaks,
  ##            colour=signal.colors[["estimate"]],linetype="dashed",lwd=1)+
  scale_fill_manual(values=breakpoint.colors)+
  ##facet_grid(segments~i,scales="free",space="free",labeller=labeller)+
  facet_grid(.~i,scales="free",space="free",labeller=labeller)+
  scale_x_continuous("position on chromosome 2 (mega base pairs)",
                     breaks=br)+
  theme_grey()
print(p)
