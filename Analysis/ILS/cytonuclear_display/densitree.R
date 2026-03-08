#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-07-22
# Description: Plot densitree to display the gene trees distinction

# ==== Main Script Start ====
library(phangorn)
library(ape)

# Load the tree
sptree <- read.tree("rosa_ags353_treeshrink_sp_rt.tre")
sptree$edge.length <- NULL

genetree <- read.tree("rosa_ags353_treeshrink_genetrees.tre")

# 清理所有树的tips，只保留物种树中的tips 
cleaned_trees <- lapply(genetree, function(tr) {
  tr$edge.length <- NULL
  tr
})

class(cleaned_trees) <- "multiPhylo"


# 如果树太多可以抽取其中一部分用于建树
# if (length(genetree) > 10) {
#   genetree_2 <- genetree[sample(1:length(genetree), 10)]
# }
# class(cleaned_trees_2) <- "multiPhylo"



### 重复物种树多次，然后读入
spmore <- replicate(2,sptree,simplify = F)
class(spmore) <- "multiPhylo"

spmore <- read.tree("./rosa_ags353_treeshrink_sp_rt_1.tre")


class(spmore) <- "multiPhylo"


# 绘制 DensiTree 图
# cladogram
# phylogram
pdf("densitree_cladogram.pdf",width=10,height=15)
densiTree(cleaned_trees,
          consensus = sptree,
          type = "cladogram",
          label.offset=2,
          width=1,
          alpha = 0.04,
          scaleX =T,
          direction = "rightwards",
          col = "#272361",
          scale.bar = F)
par(new=TRUE)
densiTree(spmore,
          consensus = sptree,
          type = "cladogram",
          label.offset=0.01,
          width=2.5,
          alpha = 1,
          scaleX =T,
          direction = "rightwards",
          col = "black",
          scale.bar = F)
dev.off()


       