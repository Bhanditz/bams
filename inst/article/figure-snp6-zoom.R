works_with_R("2.15.2",SegAnnot="1.1")

data(profiles)

pro <- profiles$hi
chr <- pro$pro
anns <- pro$ann

win <- function(min,max)data.frame(min=min*1e5,max=max*1e5)
windows <- rbind(##win(100,101),##keep these
                 win(200,201),## two dummy windows
                 win( 357.5, 359),
                 win(1061.75,1063.25))
anns$mid <- with(anns,(min+max)/2)
chr.df <- data.frame()
ann.df <- data.frame()
##showSegs <- c(13)
##seg.df <- data.frame()
##breaks <- data.frame()
for(i in 3:4){
  w <- windows[i,]
  ## sm.df <- subset(results$seg,!(first.base>w$max | last.base<w$min))
  ## sm.df$first.base[sm.df$first.base < w$min] <- w$min
  ## sm.df$last.base[sm.df$last.base > w$max] <- w$max
  ## breaks <- rbind(breaks,{
  ##   data.frame(subset(results$break.df,w$min < base & base < w$max &
  ##                     segments%in%showSegs),i)
  ## })
  ## seg.df <- rbind(seg.df,{
  ##   data.frame(subset(sm.df,segments%in%showSegs),i)
  ## })
  ann.df <- rbind(ann.df,{
    d <- subset(anns,
                (w$min < min & min < w$max)|
                (w$min < max & max < w$max))
    transform(d,i=i,
              min=ifelse(w$min>min,w$min,min),
              max=ifelse(w$max<max,w$max,max))
  })
  chr.df <- rbind(chr.df,{
    data.frame(subset(chr,w$min < position & position < w$max),i)
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
p <- ggplot()+
  geom_tallrect(aes(xmin=min/1e3,xmax=max/1e3,fill=annotation),
                colour="black",data=ann.df)+
  geom_point(aes(position/1e3,logratio),data=chr.df)+
  ## geom_segment(aes(first.base/1e3,mean,xend=last.base/1e3,yend=mean),
  ##              data=seg.df,colour=signal.colors[["estimate"]],lwd=1)+
  ## geom_vline(aes(xintercept=base/1e3),data=breaks,
  ##            colour=signal.colors[["estimate"]],linetype="dashed",lwd=1)+
  scale_fill_manual(values=breakpoint.colors)+
  facet_grid(.~i,scales="free",space="free",labeller=labeller)+
  scale_x_continuous("position on chromosome 2 (kilo base pairs)",
                     breaks=c(seq(35600,35900,by=20),seq(106130,106370,by=20)))
pdf("figure-snp6-zoom.pdf")
print(p)
dev.off()
