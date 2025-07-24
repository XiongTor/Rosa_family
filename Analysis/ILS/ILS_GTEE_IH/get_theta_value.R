#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-07-04
# Description: Get the theta value from the ASTRAL tree and the rescaled branch lengths tree.

# ====== Main Script Start ====
#R
library(ape)

args <- commandArgs(TRUE)

tree_coal <- read.tree(args[1])  # ASTRAL或MPEST的
tree_mut <- read.tree(args[2])   # IQTREE或RAxML的

# 提取边和分支长度
edge_mut <- tree_mut$edge
edge_coal <- tree_coal$edge

length_mut <- tree_mut$edge.length
length_coal <- tree_coal$edge.length

# 检查边是否对齐
if (!all(edge_mut == edge_coal)) {
  print("ERROR: The edge matrices do not match. You must fix the tree order.")
}

# 对齐边：找出每条边（子节点）在mut树中的索引
match_idx <- match(tree_mut$edge[,2], tree_coal$edge[,2])
print("Save the error")
# 注意：必须对 child（第二列）对齐，因为 tip 和内部节点编号是一致的

# 排序后的 coalescent 分支长度
length_coal_sorted <- tree_coal$edge.length[match_idx]

# θ 计算
theta <- tree_mut$edge.length / length_coal_sorted

# 把每条边的起点、终点和theta输出
df <- data.frame(
  parent = tree_mut$edge[,1],
  child = tree_mut$edge[,2],
  mutation_unit = tree_mut$edge.length,
  coalescent_unit = length_coal_sorted,
  theta = theta
)

# 可选：过滤掉 tip 边，只保留内部节点
Ntip <- length(tree_mut$tip.label)
df_internal <- df[df$child > Ntip, ]

# 写出结果
write.table(df_internal, "theta_per_node.tsv", sep="\t", quote=FALSE, row.names=FALSE)

