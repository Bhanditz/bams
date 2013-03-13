works_with_R("2.15.2",SegAnnot="1.1",ggplot2="0.9.3")

pdf("figure-snp6-both.pdf",w=9)

library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2)))

data(profiles)

pro <- profiles$hi
chr <- pro$pro

writeLines(as.character(nrow(chr)),"snp6-chr2-probes.tex")

anns <- pro$ann
win <- function(min,max)data.frame(min=min*1e5,max=max*1e5)
windows <- rbind(
                 win( 148, 171),
                 win( 355, 361),
                 win(1060,1065))
anns$mid <- with(anns,(min+max)/2)
chr.df <- data.frame()
ann.df <- data.frame()
for(i in 1:nrow(windows)){
  w <- windows[i,]
  chr.df <- rbind(chr.df,{
    data.frame(subset(chr,w$min < position & position < w$max),i)
  })
  ann.df <- rbind(ann.df,{
    data.frame(subset(anns,w$min < mid & mid < w$max),i)
  })
}
writeLines(as.character(nrow(chr.df)),"snp6-chr2-shown.tex")
labeller <- function(var,val){
  sprintf("window %d",val)
}
br <- c(6.5,7.0,seq(15,17,by=0.5),35.5,36,106,106.5)
names(br) <- as.character(br)
p <- ggplot()+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                colour="black",data=ann.df)+
  geom_point(aes(position/1e6,logratio),data=chr.df)+
  scale_fill_manual(values=breakpoint.colors)+
  facet_grid(.~i,scales="free",space="free",labeller=labeller)+
  scale_x_continuous("position on chromosome 2 (mega base pairs)",
                     breaks=br)+
  theme_grey()

pushViewport(viewport(layout.pos.row=1))
print(p,newpage=FALSE)
popViewport()

## Keep the first dummy window so that the others are numbered 2 and
## 3.
windows <- rbind(
                 win(200,201), ##dummy
                 win( 357.5, 359),
                 win(1061.75,1063.25))
anns$mid <- with(anns,(min+max)/2)
chr.df <- data.frame()
ann.df <- data.frame()
for(i in 2:3){
  w <- windows[i,]
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
p <- ggplot()+
  geom_tallrect(aes(xmin=min/1e3,xmax=max/1e3,fill=annotation),
                colour="black",data=ann.df)+
  geom_point(aes(position/1e3,logratio),data=chr.df)+
  scale_fill_manual(values=breakpoint.colors)+
  facet_grid(.~i,scales="free",space="free",labeller=labeller)+
  scale_x_continuous("position on chromosome 2 (kilo base pairs)",
                     breaks=c(seq(35600,35900,by=25),seq(106200,106400,by=25)))+
  theme_grey()
pushViewport(viewport(layout.pos.row=2))
print(p,newpage=FALSE)
popViewport()

dev.off()
