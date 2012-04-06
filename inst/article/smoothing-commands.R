data(neuroblastoma,package="neuroblastoma")
cids <- levels(neuroblastoma$profiles$profile.id)
commands <-
  sprintf("%s --vanilla --args '%s' < %s",
          R.home(file.path("bin","R")),
          cids,
       system.file(file.path("article","smooth-one.R"),package="bams"))
for(cmd in commands){
  f <- tempfile()
  qsubcmd <- sprintf("qsub %s",f)
  writeLines(cmd,f)
  system(qsubcmd)
}
