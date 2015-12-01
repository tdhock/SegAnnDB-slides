works_with_R("2.15.2",tikzDevice="0.6.3",ggplot2="0.9.3")

maketikz <- function(f,p){
  tikz(sprintf("figure-%s.tex",f),h=3,w=4.5)
  print(p)
  dev.off()
}
options(tikzMetricsDictionary="tikzMetrics")

