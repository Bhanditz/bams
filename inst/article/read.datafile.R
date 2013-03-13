read.datafile <- function
### Load a file from my server, with caching on the local filesystem.
(f,
### Name of file to download and save, relative to the current
### directory, and to the base url directory.
 base="http://cbio.ensmp.fr/~thocking/neuroblastoma"
### Directory on my server.
 ){
  tryCatch({
    read.csv(f)
  },error=function(e){
    u <- sprintf("%s/%s",base,f)
    download.file(u,f)
    read.csv(f)
  })
### The data.frame read via read.csv.
}
