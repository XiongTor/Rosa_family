#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-12-21
# Description: 

# ==== Main Script Start ====

## seqinr ------------------------------------------------------
install.packages("seqinr")
library(seqinr)
library(pheatmap)
library(grid)
library(gridExtra)
library(ggplotify)

GRP_without <- read.alignment("new_251228/PrSGRP4_protein_without.trim.fasta",format = "fasta")
GRP_all <- read.alignment("new_251228/PrSGRP_protein.trim.fasta",format = "fasta")

OPT_inter <- read.alignment("new_251228/PrSOPT_interaction.mafft",format = "fasta")
OPT_all <- read.alignment("new_251228/PrSOPT4_protein.trim_all.fasta",format = "fasta")
OPT_100 <- read.alignment("new_251228/PrSOPT4_window_101-200.fasta",format = "fasta")
OPT_500 <- read.alignment("new_251228/PrSOPT4_window_501-600.fasta",format = "fasta")

# 开始绘图
## 设置参数
GRP <- GRP_all
OPT <- OPT_inter

## 绘图
GRP_dist <- dist.alignment(GRP, matrix = "identity" )
GRP_mat <- as.matrix(GRP_dist)
GRP_mat_upper <- GRP_mat
GRP_mat_upper[lower.tri(GRP_mat_upper)] <- NA


OPT_dist <- dist.alignment(OPT, matrix = "identity" )
OPT_mat <- as.matrix(OPT_dist)
OPT_mat_upper <- OPT_mat
OPT_mat_upper[lower.tri(OPT_mat_upper)] <- NA

p1 <- as.ggplot(pheatmap(
  GRP_mat_upper,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
  legend = TRUE,
  main = "GRP Pairwise Distance",
  na_col = "white"
))

p2 <- as.ggplot(pheatmap(
  OPT_mat_upper,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
  legend = TRUE,
  main = "OPT Pairwise Distance",
  na_col = "white"
))

# 创建第三个图 (scatter plot)
d1 <- as.vector(dist.alignment(GRP, matrix = "identity"))
d2 <- as.vector(dist.alignment(OPT, matrix = "identity"))
ct <- cor.test(d1, d2, method = "pearson")

# 将 base plot 转换为 grob
p3 <- as.ggplot(~{
  plot(d1, d2, pch = 16, col = rgb(0, 0, 0, 0.5),
       xlab = "GRP pairwise distance",
       ylab = "OPT pairwise distance",
       main = "MirrorTree Correlation")
  abline(lm(d2 ~ d1), col = "red", lwd = 2)
  legend("topleft",
         legend = c(paste("r =", round(ct$estimate, 3)),
                    paste("p =", signif(ct$p.value, 3))),
         bty = "n")
})


pdf("new_251228/GRP_OPT_inter.pdf",width = 15,height = 5)
  grid.arrange(p1, p2, p3, ncol = 3)
dev.off()