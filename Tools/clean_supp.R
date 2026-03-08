#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-07-03
# Description: 

# ==== Main Script Start ====

#!/usr/bin/env Rscript

# 加载 ape 包
library(ape)

# 获取命令行参数
args <- commandArgs(TRUE)

infile <- args[1]
outfile <- args[2]

# 读入树，移除内部节点标签
tree <- read.tree(infile)
tree$node.label <- NULL

# 写出新的树文件
write.tree(tree, file = outfile)