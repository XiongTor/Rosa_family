# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.16
# This script is used to check the monophyly of each clade

file<-commandArgs(TRUE)

library(ape)
library(dplyr)
library(phytools)

if(!dir.exists("mono")){
  dir.create("mono")
}

#read tree
tree <- read.tree(file)
tip <- as.data.frame(tree$tip.label)
colnames(tip) <- c("tiplables")

#get the subfamily
tip$subf <- sapply(strsplit(as.character(tip$tiplables), "_"),`[`,1)

#get the subfamily name list
subfamily <- sapply(strsplit(as.character(tip$tiplables), "_"),`[`,1) %>% unique()

df <- data.frame()

for(name in subfamily){
  tips <- tip$tiplables[tip$subf == name]
  mono <- is.monophyletic(tree,tips,reroot = !is.rooted(tree))
  if (mono == FALSE && length(tips)>2 ){
    print(paste0(name,"is a not a monophyletic clade"))
    pdf(paste("./mono/",name,".pdf",sep = ""),height = 100,width = 100)
    is.monophyletic(tree, tips,reroot = !is.rooted(tree),plot=T,cex=10,lwd=50)
    dev.off()
    write.table(tips,paste("./mono/",name,".txt",sep = ""),row.names = F,sep = "\t",quote =F)
  }
}