read.datafile <- function(x){
  u <- sprintf("http://cbio.ensmp.fr/~thocking/neuroblastoma/%s",x)
  if(!file.exists(x)){
    download.file(u,x)
  }
  read.csv(x)
}
new.anns <- read.datafile("annotations.csv")
names(new.anns)[names(new.anns)=="profile_id"] <- "profile.id"
all.new <- subset(new.anns,type=="breakpoints")
data(neuroblastoma,package="neuroblastoma")
neuroblastomaDetailed <- with(all.new,{
  d <- data.frame(profile.id=factor(profile.id,levels(neuroblastoma$pro$pro)),
                  chromosome=factor(chromosome,c(1:22,"X","Y")),
                  min,max,annotation=as.character(annotation))
  subset(d, !is.na(profile.id))
})
nas <- sapply(neuroblastomaDetailed,function(x)sum(is.na(x)))
stopifnot(all(nas==0))
save(neuroblastomaDetailed,file="../../data/neuroblastomaDetailed.RData")
