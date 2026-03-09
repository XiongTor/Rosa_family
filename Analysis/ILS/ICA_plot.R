#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2026-03-09
# Description: plot IC/ICA result

# ==== Main Script Start ====


library(stringr)
library(ape)
library(ggtree)
library(cowplot)
library(ggplot2)

file <- "RAxML_Corrected_Probabilistic_IC_Score_BranchLabels.T4"

raw <- readLines(file, warn = FALSE)

# 1.调整树结构，使得ICA和IC值能放到支持率的位置，使其能再系统发育树中显示出来   ##############
# 捕获分支长度、IC、ICA
pattern <- ":([0-9eE\\.\\-]+)\\[([^,]+),([^\\]]+)\\]"
replacement_ic <- "\\2:\\1 "
replacement_ica <- "\\3:\\1 "

cleaned_ic <- str_replace_all(raw, pattern, replacement_ic)
cleaned_ica <- str_replace_all(raw, pattern, replacement_ica)

writeLines(cleaned_ic, "tree_ic_cleaned.nwk")
writeLines(cleaned_ica, "tree_ica_cleaned.nwk")

tree_ic <- read.tree("tree_ic_cleaned.nwk")
tree_ica <- read.tree("tree_ica_cleaned.nwk")



# 2. 画图展示：
# 2.1 封装绘图函数 ==========
plot_ic_tree <- function(tree, title, t = 0.5) {
  
  # 提取节点标签（支持值）
  node_support <- as.numeric(tree$node.label)
  
  # 创建数据框
  df <- data.frame(
    node = (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree)),
    support = node_support,
    stringsAsFactors = FALSE
  )
  
  # NA值设为NA（不显示）而不是0
  df$support[is.na(df$support)] <- NA
  
  # 关键修复：先创建基础树，再合并数据，避免fortify警告
  p <- ggtree(tree)
  p <- p %<+% df
  
  # 使用ifelse替代cut，避免产生额外的属性
  support_vec <- p$data$support
  
  p$data$group <- ifelse(
    is.na(support_vec), 
    "NA",
    ifelse(
      support_vec < -t, 
      paste0("Conflicting (< ", -t, ")"),
      ifelse(
        support_vec < 0,
        paste0("Weak [", -t, ", 0)"),
        ifelse(
          support_vec <= t,
          paste0("Moderate [0, ", t, "]"),
          paste0("Strong (> ", t, ")")
        )
      )
    )
  )
  
  # 转换为因子，固定顺序
  group_levels <- c(
    paste0("Conflicting (< ", -t, ")"),
    paste0("Weak [", -t, ", 0)"),
    paste0("Moderate [0, ", t, "]"),
    paste0("Strong (> ", t, ")"),
    "NA"
  )
  p$data$group <- factor(p$data$group, levels = group_levels)
  
  # 定义颜色和大小
  cols <- c("red", "orange", "steelblue", "darkblue", "gray80")
  sizes <- c(5, 2.5, 2.5, 5, 1)
  names(cols) <- group_levels
  names(sizes) <- group_levels
  
  # 关键修复：使用geom_nodepoint替代geom_point2，更稳定
  p_plot <- p +
    geom_nodepoint(
      aes(
        size = group,
        color = group
      ),
      shape = 16,
      alpha = 0.7,
      na.rm = TRUE  # 自动移除NA值，避免警告
    ) +
    geom_tiplab(size = 2.5) +
    scale_color_manual(
      name = "Support",
      values = cols,
      breaks = group_levels[group_levels != "NA"]  # 不显示NA在图例中
    ) +
    scale_size_manual(
      name = "Support",
      values = sizes,
      breaks = group_levels[group_levels != "NA"]
    ) +
    ggtitle(title) +
    theme(
      legend.position = c(0.15, 0.85),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  return(p_plot)
}

# 3. 绘制两棵树
p_plot_IC <- plot_ic_tree(tree_ic, "IC (Internode Certainty)", t = 0.5)
p_plot_ICA <- plot_ic_tree(tree_ica, "ICA (Internode Certainty All)", t = 0.5)

# 4. 并排显示
combined_plot <- plot_grid(
  p_plot_IC, 
  p_plot_ICA, 
  ncol = 2, 
  align = "h",
  labels = c("A", "B"),
  label_size = 14
)

print(combined_plot)

pdf("IC_ICA_trees.pdf", width = 12, height = 6)
  print(combined_plot)
dev.off()

