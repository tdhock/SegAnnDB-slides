works_with_R("2.15.2",bams="1.2",neuroblastoma="1.0")
if(getOption("repos")=="@CRAN@"){
  options(repos=c("http://mirror.ibcp.fr/pub/CRAN/",
            "http://cran.r-project.org"))
}
data(neuroblastoma)
clinical <- read.csv("clinical-limited.csv")
pids <- clinical$profile.id
pro <- subset(neuroblastoma$profiles,profile.id%in%pids & chromosome != "Y")
write.csv(pro,"profiles.csv",row.names=FALSE,quote=FALSE)
