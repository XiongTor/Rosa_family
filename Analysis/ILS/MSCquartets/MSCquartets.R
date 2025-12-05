# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.05.04
# This script is used to record the MSCquartets analysis
# install.packages("MSCquartets")

library(ape)
library(phangorn)
library(MSCquartets)

### Reading the tree file and couning the quartets
sptree=read.tree("rosa_ags353_sptree_rt_collapse10.tre")
sptree2 <- write.tree(sptree)
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"

gtrees=read.tree("rosa_ags353_genetrees_collapse10.tre")

### 得到数据集中所有物种名字
tnames=taxonNames(sptree_list)
   

# 统计基因树中所有可能的quartet的三种拓扑出现次数
QT=quartetTableParallel(gtrees,tnames,numCores=10)
#dim(QT)
save(QT,file = "QT.RData")

# 将QT中未分辨的计数分配给三种拓扑或者直接丢弃
RQT=quartetTableResolved(QT)
save(RQT,file = "RQT.RData")
### 根据共祖模型（multispecies coalescent, MSC）判断拓扑的分布情况是否与理论值显著不符
# 对应的α值
pTable=quartetTreeTestInd(RQT,model="T1",speciestree=sptree2)
pTable_2=quartetTreeTestInd(RQT,model="T3")

# 同上，但是是检验拓扑分布是否符合MSC star模型（1/3，1/3，1/3）预期的分布
# 对应的β值
Qtest=quartetStarTestInd(pTable)
Qtest_2=quartetStarTestInd(pTable_2)

#保存数据与再次读取数据
save(QT, RQT, pTable,pTable_2, Qtest,Qtest_2, file = "quartet_analysis_results.RData")
#load("quartet_analysis_results.RData")


### plot
pdf("MSCquartets_result_0.01_T3.pdf",width=16,height=12)
par(mfrow=c(1,2))
  quartetTestPlot(Qtest, "T3", alpha=.01, cex = 0.5,  beta=.05) 
  quartetTestPlot(Qtest, "T3", alpha=.01, cex = 0.5,  beta=.01) 
dev.off()

pdf("MSCquartets_result_0.001_T3.pdf",width=16,height=12)
par(mfrow=c(1,2))
quartetTestPlot(Qtest, "T3", alpha=.001, cex = 0.5,  beta=.05) 
quartetTestPlot(Qtest, "T3", alpha=.001, cex = 0.5,  beta=.01) 
dev.off()

pdf("MSCquartets_result_0.0001_T3.pdf",width=16,height=12)
par(mfrow=c(1,2))
quartetTestPlot(Qtest, "T3", alpha=.0001, cex = 0.5, beta=.05) 
quartetTestPlot(Qtest, "T3", alpha=.0001, cex = 0.5, beta=.01) 
dev.off()

pdf("MSCquartets_result_0.00001_T3.pdf",width=16,height=12)
par(mfrow=c(1,2))
quartetTestPlot(Qtest, "T3", alpha=.00001,  cex = 0.5, beta=.05)
quartetTestPlot(Qtest, "T3", alpha=.00001,  cex = 0.5, beta=.01) 
dev.off()


# count the 
# 阈值设定
alpha_val <- 0.00001
beta_val  <- 0.05

# 提取列
p_T1   <- Qtest[, "p_T3"]
p_star <- Qtest[, "p_star"]

# 去掉 NA
p_T1   <- p_T1[!is.na(p_T1)]
p_star <- p_star[!is.na(p_star)]

# 总数
n_total <- length(p_T1)

# 各类统计
n_alpha <- sum(p_T1 < alpha_val)
n_beta  <- sum(p_star > beta_val)
n_both  <- sum(p_T1 < alpha_val & p_star > beta_val)

# 比例
alpha_ratio <- n_alpha / n_total
beta_ratio  <- n_beta / n_total
both_ratio  <- n_both / n_total

cat(sprintf("超过 α 阈值 (p_T1 < %.5f): %d / %d = %.2f%%\n", alpha_val, n_alpha, n_total, alpha_ratio*100))
cat(sprintf("超过 β 阈值 (p_star > %.5f): %d / %d = %.2f%%\n", beta_val, n_beta, n_total, beta_ratio*100))
cat(sprintf("同时超过两者阈值: %d / %d = %.2f%%\n", n_both, n_total, both_ratio*100))


### inferences the tree
#stree=QDC(gtrees)
#write.tree(stree,"QDS_tree.tre")

#stree=WQDC(gtrees)
#write.tree(stree,"WQDS_tree.tre")

### plot the network tree
#NANUQ(pTable, outfile = "NANUQdist",alpha = 0.000003, beta = 0.95, plot = FALSE)
