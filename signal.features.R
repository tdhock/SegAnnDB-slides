load("signal.list.RData")
load("segmentation.list.RData")
### estimation of the variance using HALL approach
varDiff <- function(x, method='HALL'){
  n = length(x)
  if(method == 'HALL'){
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(x)
    mat[2, -n] = mat[2, -1]
    mat[3, -c(n-1, n)] = mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] = mat[4, -c(1, 2, 3)]	
    return(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3))
  }
}
kmax <- max(sapply(segmentation.list,function(x)length(x$J)))
chroms <- levels(signal.list[[1]]$chrom)
remove.inf <- function(x){
  x[!is.finite(x)] <- min(x[is.finite(x)],na.rm=TRUE)
  x
}
feature.vector <- function(pid.chr){
  signal <- signal.list[[pid.chr]]
  seg <- segmentation.list[[pid.chr]]

  J <- seg$J.est
  rss <- rep(0,kmax)
  rss[1:length(J)] <- J
  names(rss) <- sprintf("rss.%d",1:length(rss))
  log.rss <- remove.inf(log(rss))
  names(log.rss) <- sprintf("log.%s",names(rss))

  med.abs.diff <- median(abs(diff(signal$logratio)))

  n <- nrow(signal)

  mse <- rss/n
  names(mse) <- sprintf("mse.%d",1:length(mse))
  log.mse <- remove.inf(log(mse))
  names(log.mse) <- sprintf("log.%s",names(mse))

  hall <- sqrt(varDiff(signal$logratio))
  
  bases <- max(signal$position)-min(signal$position)
  probes.per.base <- n/bases
  bases.per.probe <- bases/n

  chr <- as.character(signal$chromosome[1])
  is.chrom <- as.numeric(chroms == chr)
  names(is.chrom) <- sprintf("chr%s",chroms)
  
  c(## variance estimates
    mad=med.abs.diff,
    log.mad=log(med.abs.diff),
    hall=hall,
    log.hall=log(hall),
    ## signal size
    emilie=log(2*log(n)+5),
    n=n,
    log.n=log(n),
    log2.n=log(log(n)),
    ## bases and probe density info
    bases=bases,
    log.bases=log(bases),
    probes.per.base=probes.per.base,
    log.probes.per.base=log(probes.per.base),
    bases.per.probe=bases.per.probe,
    log.bases.per.probe=log(bases.per.probe),
    ## model error
    rss,log.rss,mse,log.mse,
    ## chrom indicator
    is.chrom)
}
feature.vector("357.19")
feature.vector("349.21")
## ignore Y chrom since it is too small and causes errors with the
## feature calculation.
to.model <- grep("Y",names(signal.list),value=TRUE,invert=TRUE)
feature.vectors <- lapply(to.model,function(pid.chr){
  print(pid.chr)
  feature.vector(pid.chr)
})
## make sure they are all the same size
sizes <- sapply(feature.vectors,length)
stopifnot(all(sizes == sizes[1]))
## assemble the feature matrix
signal.features <- do.call(rbind,feature.vectors)
rownames(signal.features) <- to.model
## make sure there are no missing or infinite values in the features
bad.signals <- apply(!is.finite(signal.features),1,any)
print(signal.features[bad.signals,,drop=FALSE])
stopifnot(all(bad.signals == FALSE))
feature.sd <- apply(signal.features,2,sd)
signal.features <- signal.features[,feature.sd > 0]
save(signal.features,file="signal.features.RData")
