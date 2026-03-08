# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.18
# This script is used to check the monophyly of each clade by R package MonoPhy
# usage: Rscript monophy.R $treefile

file<-commandArgs(TRUE)

library(MonoPhy)

if(!dir.exists("monophy")){
  dir.create("monophy")
}

print(file[1])
tree<-read.tree(file[1])

# print(file[2])
# class<-read.csv(file[2],header = F)

#AssessMonophyly 
solution0 <- AssessMonophyly(tree)

#Summary
mono_clade_number<-GetSummaryMonophyly(solution0)

#Result
outlier<-GetOutlierTips(solution0, taxa = NULL, taxlevels='ALL')
df<-stack(outlier$Genera)

pdf("./monophy/result_plot.pdf",height = 200,width = 200)
  PlotMonophyly(solution0, tree,type="fan",plot.type='monophyly',label.offset=10,edge.width=3,ladderize=TRUE, cex=2)
dev.off()

write.csv(mono_clade_number,"monophy/mono_clade_number.csv",quote = F,row.names = FALSE)
write.csv(df,"monophy/outlier.csv",quote = F,row.names = FALSE)