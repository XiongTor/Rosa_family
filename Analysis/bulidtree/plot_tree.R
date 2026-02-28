#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2026-02-09
# Description: to auto plot tree for project rosaceae_cytonuclear

# ==== Main Script Start ====


# 加载必要的包
library(ape)
library(ggtree)
library(dplyr)
library(ggsci)
library(ggplot2)

# 1. 数据准备与参数设置 ------------------------------------------------------------------------
## 生成树与分组
tree <- read.tree("tree/rosa_orthofinder_MO_treeshrink_sp_rt_rename.tre")

## 分组信息
group_list <- read.csv("Rosaceae_genus_accepted_lastest.csv")
 
## 与树中的tips一一对应
tips <- as.data.frame(tree$tip.label)
colnames(tips) <- "species"
tips$genus <- sapply(strsplit(tips$species, "_"), `[`, 1)

groups <- merge(tips, group_list, by.x = "genus", by.y = "Genera", all.x = TRUE) %>% 
  select(genus, species, Tribes, Subfamilies)
groups$Tribes[is.na(groups$Tribes)] <- "outgroups"
groups$Subfamilies[is.na(groups$Subfamilies)] <- "outgroups"


## 提取族作为基本的分组信息
group_info <- data.frame(
  label = groups$species, 
  group = groups$Tribes,
  subfamily = groups$Subfamilies
)

## 给不同的分组定义颜色
my_colors <- setNames(
  pal_d3("category20")(length(unique(group_info$group))),
  unique(group_info$group)
)

## 给不同的亚科定义颜色（用于竖线）
subfamily_colors <- setNames(
  pal_d3("category10")(length(unique(group_info$subfamily))),
  unique(group_info$subfamily)
)


# 2. 自定义函数及正式画图 ------------------------------------------------------------------
## 获取树数据
get_tree_data <- function(tree) {
  p <- ggtree(tree, branch.length = "none")
  return(p$data) 
}

## 获取分组块信息
get_blocks <- function(data, y_col) {
  data <- data %>% arrange(!!sym(y_col))
  data$block_id <- cumsum(c(TRUE, data$group[-1] != data$group[-nrow(data)]))
  return(data)
}

## 生成渐变矩形数据
calc_gradient_rects <- function(block_df, x_limit, y_col, reverse_direction = FALSE) {
  gradient_dat <- data.frame()
  unique_blocks <- unique(block_df$block_id)
  n_layers <- 500
  
  for(bid in unique_blocks) {
    sub_df <- block_df %>% filter(block_id == bid)
    current_group <- sub_df$group[1]
    
    y_vals <- sub_df[[y_col]]
    ymin <- min(y_vals) - 0.5
    ymax <- max(y_vals) + 0.5
    
    if(reverse_direction) {
      x_start <- -x_limit * 0.3
      x_end <- x_limit
    } else {
      x_start <- -x_limit * 0.3
      x_end <- x_limit
    }
    
    x_seq <- seq(x_start, x_end, length.out = n_layers + 1)
    
    for(i in 1:n_layers) {
      progress <- i / n_layers
      alpha_val <- (progress^1.5) * 1
      
      # 关键：让每个矩形稍微重叠，消除缝隙
      overlap <- (x_seq[2] - x_seq[1]) * 0# 1%的重叠
      
      chunk <- data.frame(
        xmin = x_seq[i],
        xmax = x_seq[i + 1] + overlap,  # 向右延伸一点点
        ymin = ymin,
        ymax = ymax,
        alpha = alpha_val,
        group = current_group
      )
      
      gradient_dat <- rbind(gradient_dat, chunk)
    }
  }
  return(gradient_dat)
}

# 3. 计算数据并绘图 ----------------------------------------------------------------------------
## 获取树数据
data <- get_tree_data(tree)

## 画布边界---需要调节
max_x1 <- max(data$x, na.rm = TRUE)
tip_offset <- 13
limit_x1 <- max_x1 + tip_offset

## 合并分组信息
combined_data <- data %>%
  left_join(group_info, by = "label") %>%
  filter(!is.na(group))

## 计算分组块和渐变数据
blocks_left <- get_blocks(combined_data, "y")
grad_left <- calc_gradient_rects(blocks_left, limit_x1, "y", reverse_direction = TRUE)

## 计算每个分组的中心位置用于添加标签
group_labels <- blocks_left %>%
  group_by(group) %>%
  summarise(
    y_center = mean(y),
    .groups = 'drop'
  )

## 计算每个亚科的范围用于添加竖线和亚科标签
subfamily_ranges <- blocks_left %>%
  group_by(subfamily) %>%
  summarise(
    y_min = min(y) - 0.5,
    y_max = max(y) + 0.5,
    y_center = mean(y),
    .groups = 'drop'
  )

# 计算竖线和标签的位置
line_x <- limit_x1 + 10 # 竖线的x位置
subfamily_label_x <- line_x + 0.5  # 亚科标签的x位置（竖线右侧）

## --- 绘制左图 --- 
p1 <- ggtree(tree, branch.length = "none") + 
  # 添加渐变背景（关键：移除边框颜色，增加线宽）
  geom_rect(data = grad_left, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                fill = group, alpha = alpha), 
            color = NA,  # 移除边框，消除黑色线条
            linewidth = 0,  # 边框宽度设为0
            inherit.aes = FALSE) +
  # 添加支持率（如果树中有bootstrap值）
  geom_nodelab(aes(label = sprintf("%.2f", as.numeric(label)), 
                   subset = !is.na(as.numeric(label)) & as.numeric(label) > 0.1),
               size = 4.5, hjust = 1.2, vjust = -0.5, color = "black") +
  # 添加tip labels
  geom_tiplab(size = 3, offset = 0.5, hjust = 0) +
  # 添加分组名称标签
  geom_text(data = group_labels,
            aes(x = limit_x1 + 0.2, y = y_center, label = group),
            hjust = 0, size = 6, fontface = "bold", inherit.aes = FALSE) +
  # 添加彩色竖线（每个亚科不同颜色）
  geom_segment(data = subfamily_ranges,
               aes(x = line_x, xend = line_x, y = y_min, yend = y_max, color = subfamily),
               linewidth = 1.2, inherit.aes = FALSE) +
  scale_color_manual(values = subfamily_colors) +
  # 添加亚科名称标签
  geom_text(data = subfamily_ranges,
            aes(x = subfamily_label_x, y = y_center, label = subfamily),
            hjust = 0, size = 8, fontface = "bold", color = "gray20", inherit.aes = FALSE) +
  scale_fill_manual(values = my_colors) +
  scale_alpha_identity() + 
  coord_cartesian(xlim = c(-limit_x1 * 0.3, limit_x1 + 10), expand = FALSE, clip = "off") +
  # theme_tree2() +
  theme(
    legend.position = "none", 
    plot.margin = margin(5, 180, 5, -50),
    panel.background = element_blank()  # 确保背景干净
  )

# 显示图形
print(p1)

# 如果需要保存
# ggsave("phylo_tree_smooth.png", p1, 
#        width = 12, height = 15,
#        device = cairo_pdf,  # 使用cairo设备，渲染更平滑
#        dpi = 300)

# tiff("phylo_tree_smooth.tiff", 
#      width = 12, height = 15, units = "in", 
#      res = 300, compression = "lzw")
# print(p1)
# dev.off()

png("phylo_tree_high_res.png", width = 13, height = 15, units = "in", res = 600)
print(p1)
dev.off()
