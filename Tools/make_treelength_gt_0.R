# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.04.18
# This script is used to make the tree length greater than 0

file<-commandArgs(TRUE)

library(ape)
library(dispRity)

tree<-read.tree(file)

## Removing zero branch lengths
tree_good <- dispRity::remove.zero.brlen(tree)

write.tree(tree_good,"tree_without_zero_length.tre")
## Checking zero branch lengths
any(tree_good$edge.length == 0)