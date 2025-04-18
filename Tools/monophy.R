# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.18
# This script is used to check the monophyly of each clade by R package MonoPhy

tree<-commandArgs(TRUE)

library(MonoPhy)

if(!dir.exists("monophy")){
  dir.create("monophy")
}


phy <- read.tree(file=tree)
#AssessMonophyly 
solution0 <- AssessMonophyly(phy)
#Summary
GetSummaryMonophyly(solution0)
#Result
GetResultMonophyly(solution0)
PlotMonophyly(solution0, phy, plot.type='monophyly', monocoll=FALSE, label.offset=1, type="fan", cex=0.05, edge.width=1.5, PDF=TRUE, PDF_filename='monophy/result.pdf')