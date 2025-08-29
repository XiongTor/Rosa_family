#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-08-27
# Description: Extract ASTRAL support values with node IDs

# ==== Main Script Start ====

suppressMessages(library(ape))

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Using: Rscript extract_astral_support.R input_tree [output_csv]")
}

input_file <- args[1]
output_file <- ifelse(length(args) >= 2, args[2], "astral_support_gcf_scf.csv")

# 读取树
tree <- read.tree(input_file)

# 提取支持率（内部节点标签）
supports <- tree$node.label

# 节点编号（内部节点编号 = Ntip+1 : Ntip+Nnode）
node_ids <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)

# 组合结果
df <- data.frame(Node = node_ids, Support = supports, stringsAsFactors = FALSE)

# 去掉空值/NA
df <- df[!(is.na(df$Support) | df$Support == ""), ]

# 导出到 CSV
write.csv(df, output_file, row.names = FALSE, quote = FALSE)
