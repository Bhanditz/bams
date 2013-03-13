base <- "http://cbio.ensmp.fr/~thocking/neuroblastoma"
f <- "segmentation.list.RData"
u <- sprintf("%s/%s",base,f)
download.file(u,f)
