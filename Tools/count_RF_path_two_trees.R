#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-10-08
# Description: This script is used to count the RF and path distance between two trees

# ==== Main Script Start ====
library(ape)
library(phangorn)
library(ggplot2)

file <- commandArgs(TRUE)

sp <- read.tree(file[1])
gt <- read.tree(file[2])

rf_value <- numeric(length(gt))
path_value <- numeric(length(gt))


common <- intersect(sp$tip.label, gt$tip.label)
sp2 <- drop.tip(sp, setdiff(sp$tip.label, common))
gt2 <- drop.tip(gt, setdiff(gt$tip.label, common))
  
rf_value <- RF.dist(sp2, gt2, check.labels = TRUE,normalize = TRUE)
path_value <- path.dist(sp2, gt2, check.labels = TRUE, use.weight = FALSE)



result <- data.frame(
  RF_Distance = rf_value, 
  Path_Distance = path_value
  )

write.csv(result,"Rf_path_distance_sp_gene.csv",row.names = FALSE,quote = FALSE)
