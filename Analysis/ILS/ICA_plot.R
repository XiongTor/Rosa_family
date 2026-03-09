library(stringr)
library(ape)
library(ggtree)
library(cowplot)
library(ggplot2)

file <- "RAxML_Corrected_Probabilistic_IC_Score_BranchLabels.T4"

raw <- readLines(file, warn = FALSE)

###############  调整树结构，使得ICA和IC值能放到支持率的位置，使其能再系统发育树中显示出来   ##############
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

#计算IC/ICA
# 提取节点支持率
ic_vals  <- as.numeric(tree_ic$node.label)
ica_vals <- as.numeric(tree_ica$node.label)

# 计算 IC / ICA
ratio_vals <- ic_vals / ica_vals

# 如果某些 ICA 为0，会出现 Inf 或 NA 替换掉
ratio_vals[!is.finite(ratio_vals)] <- NA

# 把结果写回 tree_ic 或创建新树
tree_ratio <- tree_ic
tree_ratio$node.label <- ratio_vals


#检查是否新的树可以看到IC和ICA值
p <- ggtree(tree_ratio) + geom_tree() + geom_tiplab() +
  geom_nodelab(size=2.5)  # 节点支持显示 IC/ICA
p


## 画图展示：
node_support <- as.numeric(tree_ratio$node.label)

#提取出ICA,IC值
df <- data.frame(
  node = (Ntip(tree_ratio) + 1):(Ntip(tree_ratio) + Nnode(tree_ratio)),
  support = node_support,
  stringsAsFactors = FALSE
)

df$support[is.na(df$support)] <- 0

# 3) 使用 %<+% 将 df 附加到 ggtree 的数据上
p <- ggtree(tree_ratio) %<+% df

# 4) 绘图：按 support 控制点大小；只在内部节点(!isTip)画点
t <- 0.5  # 你的阈值，可改为 0.1 / 1 / 其它

# 分类：四档
p$data$group <- cut(
  p$data$support,
  breaks = c(-Inf, -t, 0, t, Inf),
  labels = c(
    paste0("< ", -t),
    paste0("[", -t, ", 0)"),
    paste0("[0, ", t, "]"),
    paste0("> ", t)
  ),
  right = TRUE
)

labels <- c(
  paste0("< ", -t),
  paste0("[", -t, ", 0)"),
  paste0("[0, ", t, "]"),
  paste0("> ", t)
)

# 对应颜色
cols <- c("red", "red", "blue", "blue")
sizes <- c(6, 3, 3, 6)

# 用 setNames
p_plot_IC <- p +
  geom_point2(
    aes(
      subset = !isTip,
      size = group,
      color = group
    ),
    shape = 16,
    fill = NA,
    alpha = 0.5
  ) +
  geom_tiplab(size = 3) +
  scale_color_manual(values = setNames(cols, labels)) +
  scale_size_manual(values = setNames(sizes, labels)) +
  # theme_tree2() +
  ggtitle("IC value")+
  theme(legend.position = c(0.1,0.9))

print(p_plot_IC)

plot_grid(p_plot_IC, p_plot_ICA, ncol = 2, align = "h")



# 5) 如需把支持值也显示为文字
p_plot + geom_text2(aes(subset = !isTip, label = round(support, 2)), hjust = -0.2, size = 3)

