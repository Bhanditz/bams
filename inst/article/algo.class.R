algo.list <-
  list("optimization-based models"=c("cghseg.mBIC","pelt.n","cghseg.k",
         "flsa.norm","flsa","pelt.default"),
       "approximate optimization"=c("dnacopy.alpha","dnacopy.default",
         "dnacopy.prune","dnacopy.sd","gada","cghFLasso","gada.default"),
       glad=c("glad.haarseg","glad.default","glad.MinBkpWeight",
         "glad.lambdabreak"))
algo.class <- c()
for(klass in names(algo.list)){
  algo.class[algo.list[[klass]]] <- klass
}
algo.class <- factor(algo.class,names(algo.list))
