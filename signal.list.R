works_with_R("2.15.2",neuroblastoma="1.0")

data(neuroblastoma)
noY <- subset(neuroblastoma$profiles,chromosome!="Y")
signal.list <- split(noY,with(noY,list(profile.id,chromosome)),drop=TRUE)
save(signal.list,file="signal.list.RData")
