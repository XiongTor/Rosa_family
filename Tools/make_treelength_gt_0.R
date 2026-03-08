#!/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.18
# This script is used to make the tree length greater than 0

# 获取命令行参数（树文件路径）
tree_file <- commandArgs(TRUE)

# 加载必要的包
library(ape)
library(dispRity)

# 读取树文件为 phylo 对象
tree <- read.tree(tree_file)

# 去除零分支长度
tree_good <- dispRity::remove.zero.brlen(tree)

# 构造输出名字  ——  在 .tre 前面加上 _nozero
out_file <- sub("\\.tre$", "_nozero.tre", tree_file)

# 保存结果
write.tree(tree_good, out_file)

# 检查是否还有零分支长度
any(tree_good$edge.length == 0)