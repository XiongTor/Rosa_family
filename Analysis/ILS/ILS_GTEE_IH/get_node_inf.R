#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-07-04
# Description: This script is used to get the node information.It was served as a helper function for the ILS analysis.

# ==== Main Script Start ====
library(ape)

args <- commandArgs(TRUE)

# 读取树
tree <- read.tree(args[1])

# 提取基本信息
Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode

# 内部节点的真实编号是：
internal_node_ids <- (Ntip + 1):(Ntip + Nnode)

# 提取对应的 node.label（顺序一一对应）
df <- data.frame(
  node = internal_node_ids,
  reticulate_index = as.numeric(tree$node.label)
)

# 写出为表格
write.table(df, file = "GTEE_err_reticulate_index.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
