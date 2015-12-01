works_with_R("2.15.2",SegAnnot="1.2",ggplot2="0.9.3.1")

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


data(profiles)
lo <- profiles$low
chr <- lo$pro
kmax <- 5
result <- with(chr,run.cghseg(logratio,position,kmax))
regions <- lo$ann
p <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  scale_fill_manual(values=breakpoint.colors)+
  geom_tallrect(aes(xmin=min/1e6,xmax=max/1e6,fill=annotation),
                data=regions,colour="grey")+
  geom_point(aes(position/1e6,logratio),data=chr, shape=1)+
  xlab("position on chromosome 17 (mega base pairs)")+
  scale_y_continuous("logratio (approximate copy number)")
height <- 4
lwd <- 1.5
pdf("figure-segannot.pdf",height=height)
print(p)
dev.off()

max.changes <-
  c("1breakpoint"=1,
    "0breakpoints"=0,
    ">0breakpoints"=Inf)
min.changes <-
  c("1breakpoint"=1,
    "0breakpoints"=0,
    ">0breakpoints"=1)
regions$max.changes <- max.changes[paste(regions$annotation)]
regions$min.changes <- min.changes[paste(regions$annotation)]

for(k in 1:kmax){
  breaks <- subset(result$break.df,segments==k)
  segs <- subset(result$segments,segments==k)
  pmodel <- p
  if(nrow(breaks)){
    pmodel <- pmodel+
      geom_vline(aes(xintercept=base/1e6),
                 linetype="dashed",colour=signal.colors["estimate"],
                 data=breaks,lwd=lwd)
  }
  error.regions.list <- list()
  for(region.i in 1:nrow(regions)){
    region <- regions[region.i, ]
    breakpoints.in.region <-
      sum(region$min < breaks$base & breaks$base < region$max)
    fp <- region$max.changes < breakpoints.in.region
    fn <- breakpoints.in.region < region$min.changes
    status <- ifelse(fp, "false positive",
                     ifelse(fn, "false negative", "correct"))
    error.regions.list[[region.i]] <-
      data.frame(region, fp, fn, status, 
                 row.names=NULL)
  }
  error.regions <- do.call(rbind, error.regions.list)
  pmodel <- pmodel+
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
    geom_segment(aes(first.base/1e6,mean,xend=last.base/1e6,yend=mean),
                 data=segs,colour=signal.colors["estimate"],lwd=lwd)
  fn <- sprintf("figure-segannot-%d.pdf",k)
  print(fn)
  pdf(fn,height=height)
  print(pmodel)
  dev.off()
}
