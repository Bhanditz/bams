path.to.R <- R.home(file.path("bin","R"))
install.cmd <- sprintf("%s CMD INSTALL ../..",path.to.R)
system(install.cmd)

data(neuroblastoma,package="neuroblastoma")
cids <- levels(neuroblastoma$profiles$profile.id)
commands <-
  sprintf("%s --vanilla --args '%s' < %s",
          path.to.R,
          cids,
       system.file(file.path("article","smooth-one.R"),package="bams"))
system(commands[1])
for(cmd in commands){
  f <- tempfile()
  qsubcmd <- sprintf("qsub %s",f)
  writeLines(cmd,f)
  system(qsubcmd)
}
