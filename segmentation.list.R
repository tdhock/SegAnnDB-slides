works_with_R("2.15.2",cghseg="1.0.1")

### Run maximum likelihood DP segmentation for the neuroblastoma
### profiles, and save the results in a form we can use later.
load("signal.list.RData")
global.kmax <- 20
cghseg.chrom <- function(chr){
  Y <- chr$logratio
  n <- length(Y)
  ##if(chr$chromosome[1]=="4")browser()
  kmax <- min(global.kmax,n)#k is the number of SEGMENTS not BREAKPOINTS
  cat(sprintf("profile=%s chrom=%s kmax=%d\n",chr$profile.id[1],
              chr$chromosome[1],kmax))
  result <- cghseg:::segmeanCO(Y,kmax)
  result$breakpoints <- lapply(Kseq <- 1:kmax,function(k){
    segment.ends <- result$t.est[k, 1:k]
    breakpoints <- segment.ends[-length(segment.ends)]
    (chr$position[breakpoints]+chr$position[breakpoints+1])/2
  })
  result$smooth <- sapply(Kseq,function(k){
    segment.ends <- result$t.est[k, 1:k]
    segment.starts <- c(1,segment.ends[-length(segment.ends)]+1)
    for(i in seq_along(segment.ends)){
      left <- segment.starts[i]
      right <- segment.ends[i]
      ##cat(sprintf("%d %d\n",left,right))
      Y[left:right] <- mean(Y[left:right])
    }
    Y
  }) ## n x kmax
  result$n <- n ## so we can calculate the model selection easier.
  result
}
one.chrom <- signal.list[[1]]
result <- cghseg.chrom(one.chrom)
plot(logratio~position,one.chrom)
lines(one.chrom$position,result$smooth[,3],col="green")
str(result)

segmentation.list <- lapply(signal.list,cghseg.chrom)

save(segmentation.list,file="segmentation.list.RData")
