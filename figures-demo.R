works_with_R("2.15.2",bams="1.7",quadmod="1.2",quadprog="1.5.4")
source("tikz.R")

### Used for plotting annotated regions.
geom_tallrect <- function(mapping=NULL, data=NULL, stat="identity", position="identity", ...){
  require(proto)
  GeomTallRect <- proto(ggplot2:::GeomRect,{
    objname <- "tallrect"
    required_aes <- c("xmin", "xmax")
    draw <- draw_groups <- function(.,data,scales,coordinates,
                                    ymin=0,ymax=1,...){
      ymin <- grid::unit(ymin,"npc")
      ymax <- grid::unit(ymax,"npc")
      with(ggplot2:::coord_transform(coordinates, data, scales),
           ggname(.$my_name(), {
        rectGrob(xmin, ymin, xmax - xmin, ymax-ymin,
                 default.units = "native", just = c("left", "bottom"), 
                 gp=gpar(
                   col=colour, fill=alpha(fill, alpha), 
                   lwd=size * .pt, lty=linetype, lineend="butt"
                   )
                 )
      }))
    }
  })
  GeomTallRect$new(mapping = mapping, data = data, stat = stat,
                   position = position, ...)
}

load("signal.list.RData")
load("signal.features.RData")
load("segmentation.list.RData")
load("exact.breakpoints.RData")
annotations <- read.csv("demo.csv")

## translate into standard notation, and filter out annotations for
## which we have no profiles
global.kmax <- max(sapply(segmentation.list,function(L)length(L$J.est)))

pid.chr <- with(annotations,sprintf("%s.%s",profile.id,chromosome))
have.profiles <- pid.chr %in% names(segmentation.list)
if(any(!have.profiles)){
  cat(set,": removing annotations with no profiles:\n")
  print(annotations[!have.profiles,])
  annotations <- annotations[have.profiles,]
}

## calculate if each annotation agrees with each model, returning
## list[[pid.chr]]matrix[annotation,k].
anns.by.pid.chr <-
  split(annotations,
        with(annotations,list(profile.id,chromosome)),
        drop=TRUE)
max.changes <-
  c("1breakpoint"=1,
    "0breakpoints"=0,
    ">0breakpoints"=Inf)
min.changes <-
  c("1breakpoint"=1,
    "0breakpoints"=0,
    ">0breakpoints"=1)
cost.matrices <- list()
errors.by.pid.chr <- list()
for(pid.chr in names(anns.by.pid.chr)){
  annotations <- anns.by.pid.chr[[pid.chr]]
  cost.mat <- matrix(NA,nrow(annotations),global.kmax)
  errors.by.ann <- list()
  for(ann.i in 1:nrow(annotations)){
    region <- annotations[ann.i,]
    pid <- as.character(region$profile.id)
    chromosome <- as.character(region$chromosome)
    results <- segmentation.list[[pid.chr]]
    if(is.null(results)){ ## profile data not available
      print(region)
      stop("profile data not available for ",pid.chr)
    }
    breakpoints.in.region <- sapply(results$breakpoints,function(positions){
      sum(region$min < positions & positions < region$max)
    })
    ann <- as.character(region$annotation)
    max.breakpoints <- max.changes[[ann]]
    min.breakpoints <- min.changes[[ann]]
    fp <- max.breakpoints < breakpoints.in.region
    fn <- breakpoints.in.region < min.breakpoints
    status <- ifelse(fp, "false positive",
                     ifelse(fn, "false negative", "correct"))
    errors.by.ann[[ann.i]] <-
      data.frame(pid.chr, region, fp, fn, status, segments=seq_along(fp),
                 row.names=NULL)
    cost.mat[ann.i,] <- fp+fn
  }
  errors.by.pid.chr[[pid.chr]] <- do.call(rbind, errors.by.ann)
  cost.matrices[[pid.chr]] <- cost.mat
}

## cost.matrices is a list of list of matrix[annotation,k] and so the
## job of this script is to translate that into a data.frame for each
## interval in the k/lambda error curve space.

bkpts <- exact.breakpoints[[1]]
cost <- cost.matrices[[1]]
exact.cost <- lapply(names(cost.matrices),function(pid.chr){
  bkpts <- exact.breakpoints[[pid.chr]]
  cost <- cost.matrices[[pid.chr]]
  exact.cost <- cost[,bkpts$K,drop=FALSE]
  intervals <- with(bkpts,{
    data.frame(min.lambda=lambda,
               max.lambda=c(lambda[-1],Inf),
               segments=K,
               cost=colSums(exact.cost))
  })
  within(intervals,{
    max.L <- log(max.lambda)
    min.L <- log(min.lambda)
    size <- max.L-min.L
  })
})
names(exact.cost) <- names(cost.matrices)

largest.run.min <- function
### Find the run of minimum cost with the largest size.
(cost,
 size
 ){
  m <- min(cost)
  is.min <- cost == m
  d <- c(diff(c(FALSE,is.min,FALSE)))
  ##print(data.frame(cost=c(cost,NA),size=c(size,NA),diff=d))
  starts <- which(d==1)
  ends <- which(d==-1)-1
  runs <- data.frame(starts,ends)
  ##print(runs)
  runs$size <- sapply(seq_along(starts),function(i){
    sum(size[ starts[i]:ends[i] ])
  })
  ##print(runs)
  largest <- which.max(runs$size)
  list(start=starts[largest],end=ends[largest])
}

## test intervals in the middle
intervals <- data.frame(cost=c(2,2,1,1,1,0,0,1,1,1,1,0,0,0,1))
intervals$size <- 1
largest <- largest.run.min(intervals$cost,intervals$size)
stopifnot(largest$start == 12)
stopifnot(largest$end == 14)

## test beginning
intervals$cost <- c(0,0,0,0,0,0,0,0,1,1,1,1,0,1,1)
largest <- largest.run.min(intervals$cost,intervals$size)
stopifnot(largest$start == 1)
stopifnot(largest$end == 8)

## test end
intervals$cost <- c(1,1,0,0,1,1,1,0,0,0,0,0,0,0,0)
largest <- largest.run.min(intervals$cost,intervals$size)
stopifnot(largest$start == 8)
stopifnot(largest$end == 15)

L.min.max <- do.call(rbind,lapply(names(exact.cost),function(pid.chr){
  intervals <- exact.cost[[pid.chr]]
  largest <- largest.run.min(intervals$cost,intervals$size)
  L <- cbind(min.L=intervals$min.L[ largest$start ],
             max.L=intervals$max.L[ largest$end ])
  rownames(L) <- pid.chr
  L
}))

pid.chr <- rownames(L.min.max)

intervals <- data.frame(L.min.max,pid.chr,signal.features[pid.chr,])

limits <- do.call(rbind,lapply(c("min","max"),function(limit){
  colname <- sprintf("%s.L",limit)
  data.frame(limit,pid.chr,signal.features[pid.chr,],
             L=L.min.max[,colname])
}))

pidlabs <- function(var,val){
  if(var=="pid.chr"){
    sprintf("signal $i=%s$",val)
  }else as.character(val)
}

## which 2 signals should we display?
finite <- subset(intervals,is.finite(max.L))
inf <- subset(intervals,!is.finite(max.L))
to.disp <- rbind(finite[finite$log.hall==max(finite$log.hall),],
                 inf[inf$log.hall==min(inf$log.hall),])
to.disp <- intervals[c("4.2", "8.5"),]

annotations <- read.csv("demo.csv")
ann.list <-
  split(annotations,with(annotations,list(profile.id,chromosome)),drop=TRUE)


profiles <- read.csv("profiles.csv")
if(packageVersion("ggplot2")<package_version(list(major=0,minor=9.2))){
  install.packages("ggplot2")
}
add.disp <- function(df){
  within(df,{
    pid.chr <- sprintf("%s.%s",profile.id,chromosome)
    disp <- ifelse(pid.chr%in%rownames(to.disp),1,1/3)
  })
}
toplot <- profiles
toplot$chromosome <- factor(toplot$chromosome,c(1:22,"X","Y"))
annotations$chromosome <- factor(annotations$chromosome,c(1:22,"X","Y"))
plabs <- function(var,val){
  if(var=="profile.id")sprintf("%s: %s",var,val)
  else as.character(val)
}

p <- ggplot()+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                data=add.disp(annotations),colour="grey")+
  geom_point(aes(position/1e6,logratio),pch=1,colour="black",
             data=add.disp(toplot))+
  scale_x_continuous("position on chromosome (mega base pairs)",
                     breaks=c(100,200))+
  facet_grid(profile.id~chromosome,scales="free_x",space="free",
             labeller=plabs)+
  theme_bw()+
  theme(panel.margin=grid::unit(0,"lines"))+
  theme(legend.position="top")+
  scale_fill_manual(values=c("1breakpoint"="#ff7d7d","0breakpoints"='#f6f4bf'))
png("figure-annotations.png",2000,1600,res=200)
print(p)
dev.off()

emph <- p+aes(alpha=disp)+scale_alpha_identity()
png("figure-annotations-emph.png",2000,1600,res=200)
print(emph)
dev.off()



## 2 signals.
sig.df <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  data.frame(signal.list[[pid.chr]],pid.chr)
}))
ann.df <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  data.frame(ann.list[[pid.chr]],pid.chr)
}))
kmax <- 5
seg.df <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  L <- segmentation.list[[pid.chr]]
  do.call(rbind,lapply(1:kmax,function(k){
    data.frame(signal.list[[pid.chr]],smooth=L$smooth[,k],pid.chr,k)
  }))
}))
brk.df <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  L <- segmentation.list[[pid.chr]]
  do.call(rbind,lapply(2:kmax,function(k){
    data.frame(position=L$breakpoints[[k]],pid.chr,k)
  }))
}))
## bump it up just a bit
ann.df$y <- max(sig.df$logratio)+1/4
sig.plot <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  ##theme(legend.position="top")+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                color="grey",
                data=ann.df)+
  geom_point(aes(position/1e6,logratio),
             shape=1,
             data=sig.df)+
  facet_grid(pid.chr~.,labeller=pidlabs)+
  ylab("logratio (approximate copy number)")+
  ##geom_text(aes((min+max)/2e6,y,label=annotation),data=ann.df,vjust=1)+
  scale_fill_manual(values=c("1breakpoint"="#ff7d7d",
                      "0breakpoints"='#f6f4bf'))+
  xlab("position on chromosome (mega base pairs)")

maketikz("2signals",sig.plot)
signal.colors <- c(estimate="#0adb0a",
                   latent="#0098ef")

disp.error.regions <- do.call(rbind, errors.by.pid.chr[rownames(to.disp)])
seg.plot <- function(segs){
  error.regions <- subset(disp.error.regions, segments==segs)
  sig.plot <- sig.plot+
    geom_tallrect(aes(xmin=min/1e6, xmax=max/1e6,
                      linetype=status),
                  color="black",
                  fill=NA,
                  size=1,
                  data=error.regions)+
    scale_linetype_manual("error type",
                          limits=c("correct", 
                            "false negative",
                            "false positive"
                                   ),
                          values=c(correct=0,
                            "false negative"=3,
                            "false positive"=1))+
    geom_line(aes(position/1e6,smooth),lwd=2,
              data=subset(seg.df,segs==k),
              colour=signal.colors["estimate"])
  if(segs>1){
    sig.plot <- sig.plot+
    geom_vline(aes(xintercept=position/1e6),lwd=1.5,lty="dashed",
               data=subset(brk.df,segs==k),colour=signal.colors["estimate"])
  }
  sig.plot
}
for(k in 1:kmax){
  p <- seg.plot(k)
  f <- sprintf("2signals%d",k)
  maketikz(f,p)
}

err.df.k <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  m <- cost.matrices[[pid.chr]]
  cost <- colSums(m)
  data.frame(segments=seq_along(cost),cost,pid.chr)
}))
err.k <- ggplot(err.df.k,aes(segments,cost))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  geom_point()+
  geom_line()+
  scale_y_continuous("annotation error $e_i(k)$\n", breaks=0:10)+
  scale_x_continuous("number of segments $k$",
                     breaks=c(1,5,10,20),
                     minor_breaks=NULL)+
  facet_grid(pid.chr~.,labeller=pidlabs, scales="free")
maketikz("err-k",err.k)

target.cost <- sapply(exact.cost,with,min(cost))
lim.df <- subset(limits,is.finite(L))
lim.fill <- c(min="white",max="black")
lim.latex <- c(min="$\\underline L_i$",max="$\\overline L_i$")
lim.df$latex <- lim.latex[as.character(lim.df$limit)]
disp <- lim.fill
names(disp) <- lim.latex

two.limits <- subset(lim.df,pid.chr%in%rownames(to.disp))
two.limits$cost <- target.cost[as.character(two.limits$pid.chr)]
err.df <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  data.frame(exact.cost[[pid.chr]],pid.chr)
}))
blank.df <- transform(subset(lim.df,select=-pid.chr),err=1)
errplot <- ggplot(err.df)+
  geom_blank(aes(L,err),data=blank.df)+
  geom_segment(aes(min.L,cost,xend=max.L,yend=cost),lwd=2)+
  ylab("annotation error $E_i(L)$\n")+
  xlab("penalty exponent $L$")+
  facet_grid(pid.chr~.,labeller=pidlabs)+
  geom_point(aes(L,cost,fill=latex),pch=21,data=two.limits)+
  scale_fill_manual("limit",values=disp)
maketikz("err",errplot)

kstar.df <- do.call(rbind,lapply(rownames(to.disp),function(pid.chr){
  d <- data.frame(exact.cost[[pid.chr]],pid.chr)
  rbind(
        transform(d,value=segments,what="segments $z^*_i(L)$"),
        transform(d,value=cost,what="error $E_i(L)$"))
}))
kstarplot <- ggplot()+
  geom_blank(aes(L,err),data=blank.df)+
  geom_segment(aes(min.L,value,xend=max.L,yend=value),lwd=2,data=kstar.df)+
  facet_grid(what~pid.chr,scales="free_y",labeller=pidlabs)+
  scale_fill_manual("limit",values=disp)+
  scale_y_continuous("",breaks=c(0,1,5,10,20),minor_breaks=NULL)+
  xlab("penalty exponent $L$")
maketikz("kstar",kstarplot)

two.limits$what <- kstar.df[nrow(kstar.df),"what"]
kstardots <- kstarplot+
  geom_point(aes(L,cost,fill=latex),pch=21,data=two.limits)+
  guides(fill="none")
maketikz("kstardots",kstardots)


targets <- ggplot()+
  geom_blank(aes(L,log.hall),data=lim.df)+
  ylab("variance estimate $\\log\\hat s_i$\n")+
  geom_segment(aes(min.L,log.hall,xend=max.L,yend=log.hall),data=to.disp)+
  geom_point(aes(L,log.hall,fill=latex),pch=21,data=two.limits)+
  scale_fill_manual("limit",values=disp)+
  xlab("penalty exponent $L$")
maketikz("2targets",targets)

makePlot <- function(layer=NULL){
ggplot()+
  layer+
  xlab("penalty exponent $L$")+
  ylab("variance estimate $\\log\\hat s_i$\n")+
  geom_point(aes(L,log.hall,fill=latex),pch=21,data=lim.df)+
  scale_fill_manual("limit",values=disp)
}

toPlot <- list(intervals=makePlot({
  geom_segment(aes(min.L,log.hall,xend=max.L,yend=log.hall),data=intervals)
}),limits=makePlot())

for(N in names(toPlot)){
  f <- sprintf("%s",N)
  maketikz(f,toPlot[[N]])
}


max.margin <- function
### Support vector interval regression for separable data. The idea is
### that we first normalize the feature matrix, giving normalized
### features x_i in R^p. Then we use a linear function f(x_i) = w'x_i
### + b to predict a log(lambda) that falls between all the log limits
### L_i^left and L_i^right and maximizes the margin. So the
### optimization problem is: max_{M,f} subject to, for all finite
### limits, L_i^right - f(x_i) >= M and f(x_i) - L_i^left >= M. Since
### we assume f is linear the problem becomes min_{w,b,M} -M subject
### to -M - w'x - b >= -L_i^right and -M + w'x + b >= L_i^left. We
### call M margin, w weights, b intercept.
(features,
### Matrix n x p of inputs: n signals, each with p features. We will
### scale these internally.
 limits,
### Matrix n x 2 of output lambda. Each row corresponds to the lower
### and upper bound of an interval on the log(lambda) which is optimal
### with respect to annotation error. Lower bound can be -Inf and
### upper bound can be Inf, which correspond to zero asymptotic
### cost. 
 verbose=0,
 ...
### ignored.
 ){
  ## reality checks.
  stopifnot(nrow(features)==nrow(limits))
  if(ncol(limits)!=2){
    cat("str(limits)=\n")
    str(limits)
    stop("limits should be a 2-column matrix")
  }
  stopifnot(is.matrix(features))
  
  ## check if there are any flat error curves, which have no limits.
  has.limits <- apply(is.finite(limits),1,any)
  ## we train the model on this subset.
  some.limits <- limits[has.limits,]
  some.features <- features[has.limits,,drop=FALSE]

  scaled <- scale(some.features)
  mu <- attr(scaled,"scaled:center")
  sigma <- attr(scaled,"scaled:scale")

  n <- nrow(scaled)
  p <- ncol(scaled)
  vars <- make.ids(margin=1,intercept=1,weights=p)
  constraints <- list(vars$margin*1 >= 0)
  for(i in 1:n){
    if(verbose >= 1)cat(sprintf("example constraints %5d / %5d",i,n))

    left <- some.limits[i,1]
    if(is.finite(left)){
      ivars <- with(vars,{
        intercept * 1 + sum(weights)*scaled[i,] + margin*-1
      })
      constraints <- c(constraints,list(ivars >= left))
    }

    right <- some.limits[i,2]
    if(is.finite(right)){
      ivars <- with(vars,{
        intercept * -1 + sum(weights)*scaled[i,]*-1 +margin*-1
      })
      constraints <- c(constraints,list(ivars >=  - right))
    }

    if(verbose >= 1)cat("\n")

  }
  const.info <- standard.form.constraints(constraints,vars)
  n.vars <- length(unlist(vars))
  Dvec <- rep(1e-10,n.vars)
  D <- diag(Dvec)
  d <- rep(0,n.vars)
  d[vars$margin] <- 1
  if(verbose >= 1)cat(sprintf("solving for %d variables and %d constraints... ",
              n.vars,length(constraints)))
  sol <- solve.QP(D,d,const.info$A,const.info$b0)
  if(verbose >= 1)cat("solved!\n")
  sol$mu <- mu
  sol$sigma <- sigma
  sol$scaled <- scaled
  sol$log.limits <- some.limits
  sol$weights <- sol$solution[vars$weights]
  sol$intercept <- sol$solution[vars$intercept]
  sol$margin <- sol$solution[vars$margin]
  ## this function will be applied to new data before applying the
  ## model.
  sol$normalize <- function(X){
    mu.mat <- matrix(mu,nrow(X),ncol(X),byrow=TRUE)
    s.mat <- matrix(sigma,nrow(X),ncol(X),byrow=TRUE)
    (X-mu.mat)/s.mat
  }
  sol$f <- function(x){
    sum(x*sol$weights)+sol$intercept
  }
  sol$predict <- function(X){
    stopifnot(is.matrix(X))
    X.norm <- sol$normalize(X)
    weights.mat <- matrix(sol$weights,nrow(X),ncol(X),byrow=TRUE)
    L.hat <- rowSums(X.norm * weights.mat) + sol$intercept
    L.hat
  }
  sol$L.pred <- apply(scaled,1,sol$f)
  sol$lambda.pred <- sol$predict(features)
  sol
### List of solver results. For a feature matrix X with p columns, you
### can use list$predict(X) to get model estimates of log(lambda).
}

data.set <- 
  list(features=as.matrix(intervals[,"log.hall",drop=FALSE]),
       limits=as.matrix(intervals[,c("min.L","max.L")]))

L <- data.set$lim


feature <- t(t(range(data.set$feat)))
colnames(feature) <- colnames(data.set$feat)

predict.df <- tryCatch({
  fit <- max.margin(data.set$feat,L,1)

  log.lambda <- fit$predict(feature)
  getdf <- function(margin,line){
    data.frame(line,
               L.min=log.lambda[1]+margin,L.max=log.lambda[2]+margin,
               f.min=feature[1],f.max=feature[2])
  }

  dist.to.separator <-
    fit$log.limits*matrix(c(-1,1),nrow(fit$log.limits),2,byrow=TRUE) -
      cbind(-fit$L.pred,fit$L.pred)
  margin <- min(dist.to.separator)
  support.vectors <- dist.to.separator-margin < 1e-3


  L.sv <- fit$log.limits[support.vectors]
  x.sv <- data.set$feat[row(fit$log.limits)[support.vectors]]
  L.regressed <- fit$predict(t(t(x.sv)))
  model.df <- rbind(getdf(0,"regression"),
                    getdf(margin,"limit"),
                    getdf(-margin,"limit"),
                    data.frame(line="margin",
                               L.min=L.sv,L.max=L.regressed,
                               f.min=x.sv,f.max=x.sv))

  model.df
},error=function(e){
  print("constraints are inconsistent, fitting the square loss model.")

  source("interval-regression.R")
  calcFit <- regression.funs$square
  fit <- calcFit(data.set$feat,L)

  log.lambda <- fit$predict(feature)
  data.frame(line="regression",
             L.min=log.lambda[1],L.max=log.lambda[2],
             f.min=feature[1],f.max=feature[2])
})





lt.vals <- c(margin=1,regression=2,limit=3)
p <- toPlot$limits+
  geom_segment(aes(L.min,f.min,xend=L.max,yend=f.max,linetype=line),
              data=predict.df)+
  scale_linetype_manual(values=lt.vals,guide="none")
maketikz("max-margin",p)
