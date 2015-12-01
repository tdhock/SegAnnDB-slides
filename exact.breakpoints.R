### This script saves the exact breakpoints in lambda for the cghseg.k
### model selection criterion. 

load("segmentation.list.RData")

lambda_sequence <- structure(function # Calculation of lambda
### Given a set of optimal cost from 1 segments to K segments, return a
### set of consecutive (k, lambda) with k being the optimal seg for
### for lambda k is the best up to the next lambda store in the list
(cost
### numeric vector: optimal costs.
 ){
  Kmax = length(cost)
  Kcurrent = Kmax
  Lcurrent = 0

  vK = Kmax
  vL = 0
  i <- 2
  while(Kcurrent > 1){
    smallerK <- 1:(Kcurrent-1)
    ## NOTE: we do not divide by n here, so we can directly interpret
    ## the coefficients that we get out of the regression model.
    cost.term <- (cost[Kcurrent] - cost[smallerK]) 
    lambdaTransition <-  cost.term / ( smallerK - Kcurrent)
    Kcurrent <- which.min(lambdaTransition)
    Lcurrent <- min(lambdaTransition)
    vL[i] <- Lcurrent
    vK[i] <- Kcurrent
    i <- i+1
  }
  breakpoints <- rep(NA,length(cost))
  breakpoints[vK] <- vL
  ## vL[i] stores the smallest lambda such that vK[i] segments is
  ## optimal. i
  list(K=vK, lambda=vL, breakpoints=breakpoints)
},ex=function(){
  ## First create some fake signal to segment.
  n=3*10^4
  x = rnorm(n) + rep(c(0, 1, 0), c(n/3, n/3, n/3))
  plot(x)
  ## Then segment it using cghseg.
  res = cghseg:::segmeanCO(x, K=100)
  ## Calculate the exact path of breakpoints in the optimal number of
  ## segments zstar.
  la = lambda_sequence(res$J.est)
  L <- log(la$lambda)
  exact.df <- data.frame(min.L=L,max.L=c(L[-1],Inf),segments=la$K)
  ## Solve the optimization using grid search.
  L.grid <- with(exact.df,seq(min(max.L)-1,max(min.L)+1,l=100))
  Kseq <- seq_along(res$J)
  lambda.grid <- exp(L.grid)
  kstar.grid <- sapply(lambda.grid,function(lambda){
    which.min(Kseq * lambda + res$J)
  })
  grid.df <- data.frame(log.lambda=log(lambda.grid),segments=kstar.grid)
  ## Compare the results.
  with(grid.df,plot(log.lambda,segments,pch=20,col="red",
    main=paste("Exact optimal segments curve (black lines)",
      "agrees with grid search (red points)")))
  exact.forplot <- exact.df
  exact.forplot$min.L[1] <- grid.df$log.lambda[1]-1
  exact.forplot$max.L[nrow(exact.forplot)] <-
    grid.df$log.lambda[nrow(grid.df)]+1
  with(exact.forplot,segments(min.L,segments,max.L,segments))
})

exact.breakpoints <- lapply(segmentation.list,function(seg.info){
  with(seg.info,lambda_sequence(J.est))
})
save(exact.breakpoints,file="exact.breakpoints.RData")
