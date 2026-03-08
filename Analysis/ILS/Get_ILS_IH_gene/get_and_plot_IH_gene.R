#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-11-19
# Description: This script is used to read the IH/ILS gene results and plot boxplots and barplots.

# ==== Main Script Start ====

library(ape)
library(phangorn)
library(ggplot2)
library(TreeDist)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

df <- read.csv("all_result.csv")

# 确保 genename 列存在
rownames(df) <- as.character(df[, 1])
df <- df[, -1]




# ========== 图1: 箱线图（按列顺序）==========
cat("\n生成箱线图...\n")

# 定义函数识别异常值
remove_outliers <- function(x) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lower <- q[1] - 1.5 * iqr
  upper <- q[2] + 1.5 * iqr
  x >= lower & x <= upper
}

# 准备数据用于箱线图
df_long <- df %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "variable", values_to = "value")

# 保持原始列顺序
df_long$variable <- factor(df_long$variable, levels = colnames(df))

# 过滤异常值
df_filtered <- df_long %>%
  group_by(variable) %>%
  filter(remove_outliers(value)) %>%
  ungroup()

# 绘制箱线图
p1 <- ggplot(df_filtered, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7,outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.3, height = 0),alpha=0.3,size=0.5) +
  #标注中位数的数值
  stat_summary(fun.data = function(x) {
    data.frame(y = max(x, na.rm = TRUE), label = round(median(x, na.rm = TRUE), 5))
  }, geom = "text", vjust = -0.5, size = 3.5, fontface = "bold") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  labs(x = "type", y = "value")

ggsave("boxplot_by_column.pdf", p1, width = 12, height = 6)
cat("箱线图已保存: boxplot_by_column.pdf\n")

# ============================== 图2: 柱状图系列 ==================================
# RF

# 读取树文件
sp <- read.tree("rosa_ags353_treeshrink_sp_rt_oneoutg_final.tre")
gt <- read.tree("rosa_ags353_genetrees_oneoutg_final.tre")

# 读取基因名字文件
gene_names <- readLines("genetrees.txt")

# 自定义函数统一标签
fix_labels <- function(sp, gt) {
  common <- intersect(sp$tip.label, gt$tip.label)
  sp_pruned <- keep.tip(sp, common)
  gt_pruned <- keep.tip(gt, common)
  list(sp = sp_pruned, gt = gt_pruned)
}

# 计算4个距离指标
n_genes <- length(gt)
rf_distances <- numeric(n_genes)
path_distances <- numeric(n_genes)

#计算RF与path距离
for (i in 1:n_genes) {
  tmp <- fix_labels(sp, gt[[i]])
  rf_distances[i] <- RF.dist(tmp$sp, tmp$gt, normalize=TRUE, rooted=FALSE)
  path_distances[i] <- path.dist(tmp$sp, tmp$gt, check.labels=TRUE, use.weight=FALSE)
}

# 计算 MDS (用于 MDC 和 MPD)
all_trees <- c(list(sp), gt)
dist_rf_matrix <- matrix(0, length(all_trees), length(all_trees))
for (i in 1:(length(all_trees)-1)) {
  for (j in (i+1):length(all_trees)) {
    tmp <- fix_labels(all_trees[[i]], all_trees[[j]])
    dist_rf_matrix[i,j] <- RF.dist(tmp$sp, tmp$gt, normalize=TRUE, rooted=FALSE)
    dist_rf_matrix[j,i] <- dist_rf_matrix[i,j]
  }
}

RF_mds <- cmdscale(as.dist(dist_rf_matrix), k=2)
centroid <- colMeans(RF_mds)
mdc_values <- sqrt(rowSums((RF_mds - centroid)^2))[-1]  # 去掉物种树
mpd_values <- rowMeans(dist_rf_matrix)[-1]  # 去掉物种树

# 整合结果
results <- data.frame(
  Line_Number = 1:n_genes,
  Gene_Name = gene_names,
  RF_distance = rf_distances,
  Path_distance = path_distances,
  MDC = mdc_values,
  MPD = mpd_values
)

# 找出每个指标的 Top 10
top10_rf <- results[order(-results$RF_distance), ][1:10, c("Line_Number", "Gene_Name", "RF_distance")]
top10_path <- results[order(-results$Path_distance), ][1:10, c("Line_Number", "Gene_Name", "Path_distance")]
top10_mdc <- results[order(-results$MDC), ][1:10, c("Line_Number", "Gene_Name", "MDC")]
top10_mpd <- results[order(-results$MPD), ][1:10, c("Line_Number", "Gene_Name", "MPD")]

# 保存结果
write.csv(top10_rf, "top10_RF_distance.csv", row.names=FALSE)
write.csv(top10_path, "top10_Path_distance.csv", row.names=FALSE)
write.csv(top10_mdc, "top10_MDC.csv", row.names=FALSE)
write.csv(top10_mpd, "top10_MPD.csv", row.names=FALSE)

### 提取IH_0.01最严重的10个基因
cat("\n生成柱状图系列...\n")

# 提取 d_IH_0.01 列，确保是数值型
d_IH_values <- as.numeric(df$IH_0.01)
names(d_IH_values) <- rownames(df)
# 提取IH_0.01最大的10个基因
top10_indices <- order(d_IH_values, decreasing = TRUE)[1:10]
top10_genes <- rownames(df)[top10_indices]




###################   开始画图   ###########################
#提取RF,path，mds,mpd最大的10个基因并画在图中

top10_rf_genes <- top10_rf$Gene_Name
top10_path_genes <- top10_path$Gene_Name
top10_mdc_genes <- top10_mdc$Gene_Name
top10_mpd_genes <- top10_mpd$Gene_Name


diff_type <- c("top10_genes","top10_rf_genes","top10_path_genes","top10_mdc_genes","top10_mpd_genes")

for(name in diff_type){
cat("\nd_IH_0.01 最高的10个基因:\n")
  
current_target_genes <- get(name)
print(data.frame(
  gene = top10_genes,
  d_IH_0.01 = d_IH_values[current_target_genes]
))

cat("\n基因总数:", length(current_gene_order), "\n")

# 准备所有柱状图
plot_list <- list()

for (col_name in colnames(df)) {
  # 提取当前列的值
  current_values <- df[, col_name]
  names(current_values) <- rownames(df)
  
  # 按当前列的值从大到小排序基因
  current_gene_order <- names(sort(current_values, decreasing = TRUE))
  
  # 准备数据
  plot_data <- data.frame(
    gene = current_gene_order,
    value = current_values[current_gene_order],
    is_top10 = current_gene_order %in% current_target_genes
  )
  
  # 保持当前列的基因顺序（从大到小）
  plot_data$gene <- factor(plot_data$gene, levels = current_gene_order)
  
  # 创建柱状图
  p <- ggplot(plot_data, aes(x = gene, y = value, fill = is_top10)) +
    geom_col(width = 1) +
    scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "gray70"),
                      labels = c("TRUE" = "Top 10 (PATH)", "FALSE" = "Others"),
                      name = "") +
    theme_test() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 10, face = "bold"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "top") +
    labs(title = col_name,
         x = "gene")
  
  plot_list[[col_name]] <- p
}

  # 保存所有柱状图到一个PDF文件
pdf(paste0(name,"_barplots_all_columns.pdf"), width = 14, height = 12)
for (i in seq(1, length(plot_list), by = 8)) {
  end_idx <- min(i + 7, length(plot_list))
  gridExtra::grid.arrange(grobs = plot_list[i:end_idx], ncol = 2, nrow = 4)
}
dev.off()
}





