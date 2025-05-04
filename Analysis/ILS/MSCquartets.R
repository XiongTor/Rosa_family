# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.05.04
# This script is used to record the MSCquartets analysis

install.packages("MSCquartets")

library(MSCquartets)

### Reading the tree file and couning the quartets
gtrees=read.tree("best_genetree/rosa_sp_gene.tre")

tnames=taxonNames(gtrees)   ###得到数据集中所有物种名字

QT=quartetTable(gtrees)

RQT=quartetTableResolved(QT)   ##转换所有quartet counts，丢弃unresolved的结果或者将其分到3个resolved计数（ mab|cd; mac|bd; mad|bc ）中。

#quartetTable(trees, taxonnames = NULL, epsilon = 0, random = 0)  构建qcCFs table：taxonnames可以指定想要分析的物种，例如tnames[1:6]取前6个物种；eposilon是最小可以认为是非0的分支长度；random是如果随机抽取4个分类群的子集，抽取的数目，如果random=0，则使用所有子集。

### Conduct hypothesis testing
pTable=quartetTreeTestInd(RQT,model="T1",speciestree=sptree)
#T3 means any topological, or you can use you species tree to do a constrained.
# warning: the species tree should be text format, not a tree object.

Qtest=quartetStarTestInd(pTable)


### plot
pdf("MSCquartets_result.pdf",width=10,height=10)
  quartetTestPlot(Qtest, "T1", alpha=.000003, beta=.05) 
dev.off()


### inferences the tree
stree=QDC(gtrees)
write.tree(stree,"QDS_tree.tre")

stree=WQDC(gtrees)
write.tree(stree,"WQDS_tree.tre")

### plot the network tree
NANUQ(pTable, outfile = "NANUQdist",alpha = 0.000003, beta = 0.95, plot = FALSE)
