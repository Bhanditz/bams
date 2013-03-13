load("segmentation.list.RData")

parameter.list <- lapply(segmentation.list[[1]],"[[","parameters")
parameter.list[sapply(parameter.list,length)==1] <- NA

save(parameter.list, file="parameter.list.RData")
