#!/usr/bin/env Rscript
# Author: Tao Xiong (Enhanced Version)
# Date: 2026-03-09
# Description: 计算物种树与基因树、基因树之间的 RF、Path 及 Quartet 距离，处理物种不全的情况。

# ==== 加载包 ====
suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggplot2)
  library(TreeDist)
  library(Quartet) # 专门用于四分体距离计算
})

# ==== 读取输入 ====
args <- commandArgs(TRUE)
if (length(args) < 2) {
  stop("使用方法: Rscript script.R <species_tree.tre> <gene_trees.tre>")
}

sp <- read.tree(args[1])
gt <- read.tree(args[2])

# 如果基因树文件中只有一棵树，将其转为列表格式保持一致性
if (class(gt) == "phylo") gt <- list(gt)

# 所有树合并，第一棵始终是物种树
all_trees <- c(list(sp), gt)
n_trees <- length(all_trees)

# ==== 初始化距离矩阵 ====
# 默认填充 NA 或最大差异
dist_rf <- matrix(0, n_trees, n_trees)
dist_path <- matrix(0, n_trees, n_trees)
dist_quartet <- matrix(0, n_trees, n_trees)

cat("开始计算距离矩阵，共计", n_trees, "棵树...\n")

# ==== 核心计算循环 ====
for (i in 1:(n_trees - 1)) {
  for (j in (i + 1):n_trees) {
    
    # 获取共有物种
    common <- intersect(all_trees[[i]]$tip.label, all_trees[[j]]$tip.label)
    n_common <- length(common)
    
    # 只有共有物种 >= 4 时，无根树的拓扑比较才有意义
    if (n_common >= 4) {
      # 剪枝
      tmp_i <- keep.tip(all_trees[[i]], common)
      tmp_j <- keep.tip(all_trees[[j]], common)
      
      # 1. 标准化 RF 距离 (处理报错的关键：显式 unroot)
      # 即使 unroot 失败也捕获错误，防止脚本崩溃
      tryCatch({
        dist_rf[i,j] <- RF.dist(unroot(tmp_i), unroot(tmp_j), normalize = TRUE)
      }, error = function(e) { dist_rf[i,j] <<- 1 })
      
      # 2. Path 距离 (进行简单的物种数标准化)
      dist_path[i,j] <- path.dist(tmp_i, tmp_j) / n_common
      
      # 3. Quartet 距离 (推荐：处理缺失物种最稳健)
      # QuartetDivergence 返回 0-1 之间的差异度
      dist_quartet[i,j] <- QuartetDivergence(QuartetStatus(tmp_i, tmp_j), similarity = FALSE)
      
    } else {
      # 共有物种太少，无法比较，赋予最大差异
      dist_rf[i,j] <- 1
      dist_path[i,j] <- NA 
      dist_quartet[i,j] <- 1
    }
    
    # 对称填充
    dist_rf[j,i] <- dist_rf[i,j]
    dist_path[j,i] <- dist_path[i,j]
    dist_quartet[j,i] <- dist_quartet[i,j]
  }
}

# 填充 Path 距离中的缺失值
dist_path[is.na(dist_path)] <- max(dist_path, na.rm = TRUE)

# ==== MDS 降维分析 ====
# 使用 add=TRUE 处理可能的非欧几里得距离问题
cat("正在进行 MDS 降维...\n")
run_mds <- function(d_matrix) {
  cmdscale(as.dist(d_matrix), k = 2, add = TRUE)$points
}

RF_ret <- run_mds(dist_rf)
Path_ret <- run_mds(dist_path)
Quartet_ret <- run_mds(dist_quartet)

# ==== 统计分散度 (MDC/MPD) ====
calc_dispersion <- function(coords, remove_outliers = TRUE) {
  centroid <- colMeans(coords)
  dist_to_centroid <- sqrt(rowSums((coords - centroid)^2))
  
  if (remove_outliers) {
    med <- median(dist_to_centroid)
    mad_val <- mad(dist_to_centroid, constant = 1)
    keep <- abs(dist_to_centroid - med) <= (3 * mad_val)
    if (sum(keep) > 2) coords <- coords[keep, ]
  }
  
  list(
    MDC = mean(sqrt(rowSums((coords - colMeans(coords))^2))),
    MPD = mean(as.matrix(dist(coords)))
  )
}

res_rf <- calc_dispersion(RF_ret)
res_path <- calc_dispersion(Path_ret)
res_quartet <- calc_dispersion(Quartet_ret)

# 保存结果
dir.create("MDS", showWarnings = FALSE)
total_res <- data.frame(
  Metric = c("RF", "Path", "Quartet"),
  MDC = c(res_rf$MDC, res_path$MDC, res_quartet$MDC),
  MPD = c(res_rf$MPD, res_path$MPD, res_quartet$MPD)
)
write.csv(total_res, "MDS/dispersion_metrics.csv", row.names = FALSE)

# ==== 绘图函数 ====
plot_mds <- function(coords, title, stats) {
  df <- data.frame(
    X = coords[,1], Y = coords[,2],
    Type = c("Species tree", rep("Gene tree", n_trees - 1))
  )
  
  label_text <- paste0("MDC = ", round(stats$MDC, 3), "\nMPD = ", round(stats$MPD, 3))
  
  ggplot(df, aes(X, Y, color = Type, shape = Type)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("Species tree" = "red", "Gene tree" = "blue")) +
    scale_shape_manual(values = c("Species tree" = 17, "Gene tree" = 1)) +
    theme_bw() +
    labs(title = title, x = "Dimension 1", y = "Dimension 2") +
    annotate("text", x = Inf, y = Inf, label = label_text, 
             hjust = 1.1, vjust = 1.5, size = 3.5, fontface = "italic") +
    theme(legend.position = "bottom")
}

# 生成三张图
p_rf <- plot_mds(RF_ret, "MDS (Normalized RF Distance)", res_rf)
p_path <- plot_mds(Path_ret, "MDS (Standardized Path Distance)", res_path)
p_quartet <- plot_mds(Quartet_ret, "MDS (Quartet Distance)", res_quartet)

ggsave("MDS/RF_MDS.pdf", p_rf, width = 8, height = 7)
ggsave("MDS/PATH_MDS.pdf", p_path, width = 8, height = 7)
ggsave("MDS/Quartet_MDS.pdf", p_quartet, width = 8, height = 7)

cat("分析完成！结果已保存至 MDS/ 目录下。\n")