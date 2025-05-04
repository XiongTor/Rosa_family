# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.18
# This script is used to check the monophyly of each clade by R package MonoPhy

file<-commandArgs(TRUE)

library(MonoPhy)

if(!dir.exists("monophy")){
  dir.create("monophy")
}

tree<-read.tree(file[1])

class<-read.table(file[2],header = F)

phy <- read.tree(file=tree)
#AssessMonophyly 
solution0 <- AssessMonophyly(tree,class)

#Summary
mono_clade_number<-GetSummaryMonophyly(solution0)

#Result
outlier<-GetOutlierTips(solution0, taxa = NULL, taxlevels='ALL')

#Intruders
intruders<-GetIntruderTips(solution0)

pdf(paste("./monophy/",name,".pdf",sep = ""),height = 10000,width = 100)
  PlotMonophyly(solution0, tree, plot.type='monophyly', ladderize=TRUE, cex=0.5)
dev.off()

write.csv(mono_clade_number,"monophy/mono_clade_number.csv",quote = F,row.names = FALSE)
write.csv(outlier,"monophy/outlier.csv",quote = F,row.names = FALSE)