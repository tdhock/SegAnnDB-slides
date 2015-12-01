works_with_R("2.15.2",neuroblastoma="1.0",ggplot2="0.9.3")
data(neuroblastoma)
clinical <- read.csv("clinical-limited.csv")
clinical$display.id <-
  sprintf("%d %s",clinical$profile.id,clinical$relapse)
clinical <- clinical[order(clinical$relapse,clinical$profile.id),]
clinical$display.id <- factor(clinical$display.id,clinical$display.id)
profiles <- subset(neuroblastoma$profiles,profile.id%in%clinical$profile.id)
merged <- merge(clinical[,c("profile.id","relapse","display.id")],profiles)
toplot <- subset(merged,chromosome!="Y")
toplot$chromosome <- factor(toplot$chromosome,c(1:22,"X","Y"))
library(grid)
p <- ggplot()+
  geom_point(aes(position/1e6,logratio),pch=1,colour="black",data=toplot)+
  scale_x_continuous("position on chromosome (mega base pairs)",
                     breaks=c(100,200))+
  facet_grid(display.id~chromosome,scales="free_x",space="free")+
  theme_bw()+
  theme(panel.margin=unit(0,"lines"))
png("figure-profiles.png",2000,1600,res=200)
print(p)
dev.off()
