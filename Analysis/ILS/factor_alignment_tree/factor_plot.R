# !/usr/bin/Rscript
# Author: XiongTao
# date: 2025.12.11
# usage: Anlayse the different factor between Chloroplast and Nuclear in Rosaceae
library(ape)
library(phangorn)
library(dplyr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(cowplot)
library(clue) 
library(ggplot2)
library(plotly)
#相关的因素，使用的因子  -------------------------------------------------------------------
# No_of_taxa
# Alignment_length
# GC_content
# Proportion_variable_sites
# Proportion_parsimony_informative
# Mean_evolutionary_rate
# Mean_evolutionary_rate_without0
# RCV
# 
# ####与系统发育树相关
# Mean_bipartition
# Mean_internal_branch
# Mean_teminal_branch
# treeness
# 
# ####与序列和树都相关
# saturation
# absolute_saturation
# treeness.RCV

# 1.数据读取 --------------------------------------------------------------------------------

## 准备叶绿体的数据
file_c <- list.files("./chloroplast/")

df_chloroplast <- list()
for (name in file_c) {
  df_chloroplast[[name]] <- read.table(paste0("./chloroplast/",name),header = T)
  df_chloroplast[[name]]$Alignment_name <- sapply(
    strsplit(df_chloroplast[[name]]$Alignment_name,"\\."),`[`,1
  )
}

# 将list中所有的数据框按照Alignment_name列进行合并
merged_chloroplast <- reduce(df_chloroplast, full_join, by="Alignment_name") %>% 
  select(Alignment_name,
         No_of_taxa,
         Alignment_length,
         GC_content,
         Proportion_variable_sites,
         Proportion_parsimony_informative,
         Mean_evolutionary_rate,
         Mean_evolutionary_rate_without0,
         RCV,
         Mean_bipartition,
         Mean_internal_branch,
         Mean_teminal_branch,
         treeness,
         saturation,
         absolute_saturation,
         treeness.RCV)

## 准备核基因的数据
file_ags <- list.files("./orthofinder/")

df_ags <- list()
for (name in file_ags) {
  df_ags[[name]] <- read.table(paste0("./orthofinder/",name),header = T)
  # df_ags[[name]]$Alignment_name <- sapply(
  #   strsplit(df_ags[[name]]$Alignment_name,"\\_"),`[`,1
  # )
}


merged_orthofinder <- reduce(df_ags, full_join, by="Alignment_name") %>% 
  select(Alignment_name,
         No_of_taxa,
         Alignment_length,
         GC_content,
         Proportion_variable_sites,
         Proportion_parsimony_informative,
         Mean_evolutionary_rate,
         Mean_evolutionary_rate_without0,
         RCV,
         Mean_bipartition,
         Mean_internal_branch,
         Mean_teminal_branch,
         treeness,
         saturation,
         absolute_saturation,
         treeness.RCV)


## 去除无用的元素
rm(file_c,file_ags,df_ags,df_chloroplast)


# 2.画图展示所有因子在叶绿体和核基因之间的数值差异 -------------------------------------------

## 先给数据框加标签方便区分来源
merged_orthofinder$dataset <- "orthofinder"
merged_chloroplast$dataset <- "chloroplast"

## 合并
df <- bind_rows(merged_orthofinder, merged_chloroplast)

## 指定要转化的列
cols <- c(
  "No_of_taxa","Alignment_length","GC_content",
  "Proportion_variable_sites","Proportion_parsimony_informative",
  "Mean_evolutionary_rate","Mean_evolutionary_rate_without0",
  "RCV","Mean_bipartition","Mean_internal_branch",
  "Mean_teminal_branch","treeness","saturation",
  "absolute_saturation","treeness.RCV"
)

## 创建一个空列表来保存每个变量的图
plot_list <- list()
n_col <- 3 #准备按照3×3的布局来进行绘图
col_first_idx <- seq(1, length(cols), by = n_col)
##循环每一列
for (i in seq_along(cols)) {
  col_name <- cols[i]
  
  # 全局 x 轴范围
  all_values <- c(merged_orthofinder[[col_name]], merged_chloroplast[[col_name]])
  x_min <- min(all_values, na.rm = TRUE)
  x_max <- max(all_values, na.rm = TRUE)
  
  p <- ggplot() +
    geom_density(data = merged_orthofinder, aes_string(x = col_name),
                 adjust = 3, fill="#2073a5", alpha=0.4) +
    geom_density(data = merged_chloroplast, aes_string(x = col_name),
                 adjust = 3, fill="#a64d6d", alpha=0.4) +
    labs(x = col_name, y = "Density") +
    xlim(x_min, x_max) +
    theme_bw()+
    theme(
      axis.title.x = element_text(size = 20),   # x 轴标题大小
      axis.title.y = element_text(size = 20),   # y 轴标题大小
      axis.text.x  = element_text(size = 18),    # x 轴刻度大小
      axis.text.y  = element_text(size = 18)     # y 轴刻度大小
    )
  
  # 判断是否需要隐藏 y 轴
  if (i == 1) {
    # 计算左上角位置
    dens_max <- max(density(merged_orthofinder[[col_name]], na.rm=TRUE)$y,
                    density(merged_chloroplast[[col_name]], na.rm=TRUE)$y)
    p <- p +
      annotate("rect", xmin = x_min, xmax = x_min + (x_max-x_min)*0.08,
               ymin = dens_max*0.85, ymax = dens_max*0.95, fill="#2073a5") +
      annotate("rect", xmin = x_min, xmax = x_min + (x_max-x_min)*0.08,
               ymin = dens_max*0.65, ymax = dens_max*0.75, fill="#a64d6d") +
      annotate("text", x = x_min + (x_max-x_min)*0.1,
               y = dens_max*0.9, label = "orthofinder", hjust = 0, size = 5) +
      annotate("text", x = x_min + (x_max-x_min)*0.1,
               y = dens_max*0.7, label = "Chloroplast", hjust = 0, size = 5)
  } else {
    # 其他图隐藏 y 轴刻度和标题
    p <- p + theme(
      axis.title.y = element_blank(),
    )
  }
  
  plot_list[[col_name]] <- p
}


## 用cowplot将所有图拼接在一起，按列排列
final <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
ggsave("all_density_plots.pdf", final, width = 15, height = 15)


## 去除无用的元素
rm(df,p,plot_list)

# 3.相关性分析，判断哪些因素相关性高 -------------------------------------------------------------

## 根据上一步的结果，挑选出叶绿体和核基因有差异的因子
diff <- c("GC_content",
  "Proportion_variable_sites","Proportion_parsimony_informative",
  "Mean_evolutionary_rate","RCV")

## 计算相关性
df <- merged_orthofinder[, cols]

res <- rcorr(as.matrix(df))

cor_mat  <- res$r   # 相关性
p_mat    <- res$P   # p value

pdf("corrplot_diff_all_factor.pdf",width = 10,height=10)
  corrplot(
    cor_mat,
    method = "square",
    tl.col = "black",
    p.mat   = p_mat,
    insig = "blank",
    addCoef.col = "black",
    number.cex = 1
  )
dev.off()

## 根据相关性的结果，将不同相关性程度的因子的数值画成直方图来表达
## 以"Proportion_variable_sites","Proportion_parsimony_informative"为例
df <- merged_orthofinder %>% 
  arrange(desc(Proportion_variable_sites)) %>%
  mutate(order=row_number()) %>% 
  select(order,all_of(diff))

df_long <- df %>%
  pivot_longer(cols = diff,
               names_to = "variable",
               values_to = "value") %>% 
  mutate(variable = factor(variable,
                           levels = c("Proportion_variable_sites",
                                      "Proportion_parsimony_informative",
                                      "Mean_evolutionary_rate",
                                      "GC_content"
                                      )))

ggplot(df_long) +
  geom_col(aes(x = order, y = value), fill="steelblue") +
  facet_grid(variable ~ ., scales = "free_y") +
  theme_bw()


col <- ggplot(df, aes(x = order)) +
  geom_col(aes(y = Proportion_variable_sites, fill="PVS"), alpha=0.5) +
  geom_col(aes(y = Proportion_parsimony_informative, fill="PPI"), alpha=0.6) +
  geom_col(aes(y = Mean_evolutionary_rate, fill="Rate"), alpha=1) +
  geom_col(aes(y = GC_content, fill="GC"), alpha=0.6) +
  scale_fill_manual(
    breaks = c("PVS", "PPI", "Rate", "GC"),   # <- 明确指定图例顺序
    values = c(
      "PVS"="#0072B2",
      "PPI"="red",
      "Rate"="black",
      "GC"="#DAA520"
    ),
    labels = c(
      "PVS"="Variable sites",
      "PPI"="Parsimony informative",
      "Rate"="Evolutionary rate",
      "GC"="GC content"
    )
  ) +
  labs(x="Genes", y="Value") +
  theme_bw() +
  theme(
    legend.position = c(0.88,0.85),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_blank()
  )

ggsave("col.pdf", col, width = 7, height = 5)


# 3.2.PCA -------------------------------------------------------------------------------------------------------------

## 添加基因类型标签
merged_orthofinder$type <- "nuclear"
merged_chloroplast$type <- "chloroplast"

## 合并
df_all <- bind_rows(merged_orthofinder, merged_chloroplast)

## 查看有哪些属性列（排除 gene_name 和 type）
setdiff(colnames(df_all), c("Alignment_name", "type"))

## 选属性列（你按实际列名修改）
props <- c("GC_content",
           "Proportion_variable_sites","Proportion_parsimony_informative",
           "Mean_evolutionary_rate","RCV")

## PCA 输入矩阵
X <- scale(df_all[, props])

## PCA计算，计算每个轴中各个因素的占比贡献
pca <- prcomp(X, center = TRUE, scale. = TRUE)
saveRDS(df_all, file = "df_all.rds")
saveRDS(pca,    file = "pca.rds")

## 统计每一个轴的解释度
summary(pca)$importance[2, ]

pca_df <- data.frame(
  Alignment_name = df_all$Alignment_name,
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  type = df_all$type
)

## 提取 PCA 映射表
pc_map <- pca_df[, c("Alignment_name", "PC1", "PC2", "PC3")]

## 按 Alignment_name merge 回核基因
merged_orthofinder <- merge(
  merged_orthofinder,
  pc_map,
  by = "Alignment_name",
  all.x = TRUE
)

## 按 Alignment_name merge 回质体基因
merged_chloroplast <- merge(
  merged_chloroplast,
  pc_map,
  by = "Alignment_name",
  all.x = TRUE
)

# 按照二维画出不同PC轴组成的平面的点的分布
p3 <- ggplot(pca_df, aes(PC2, PC3, color = type)) + 
  geom_point(size = 3, alpha=0.8) +
  stat_ellipse(level = 0.95, size=1) +
  theme_bw(base_size = 16)+
  theme(legend.position = c(0.15,0.9),
        legend.background = element_blank())


pdf("pca.pdf",width = 18,height=6)
plot_grid(p1,p2,p3,ncol=3)
dev.off()


## 三维

pca_df_3d <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  type = df_all$type
)

plot_ly(
  data = pca_df_3d,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~type,
  colors = c("#1f77b4", "#ff7f0e"),   # 两组颜色，可自行改
  marker = list(size = 4)
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

# 4.针对上述三个因素挑选高中低三个指标的基因，并计算高中低三个区间内基因树彼此之间的RF距离 ---------------------------------------------------------

# 自定义函数计算RF距离

## 自定义函数-------------
pairwise_rf_onecol <- function(
    genes,
    all_tree,
    min_tips = 4,
    pattern_suffix = "_",
    normalize = TRUE,
    rooted = FALSE
) {
  require(ape)
  require(phangorn)
  
  # 2. 匹配基因对应的树文件
  matched_files <- all_tree[
    sapply(all_tree, function(f)
      any(startsWith(basename(f), paste0(genes, pattern_suffix)))
    )
  ]
  
  if (length(matched_files) < 2) {
    warning("Matched tree files < 2, RF not computed.")
    return(matrix(NA, ncol = 1, dimnames = list(NULL, "RF")))
  }
  
  # 3. 读入树
  gene_trees <- lapply(matched_files, read.tree)
  names(gene_trees) <- sub("_.*$", "", basename(matched_files))
  
  n <- length(gene_trees)
  
  # 4. 计算 pairwise RF
  rf_vec <- c()
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      common_tips <- intersect(
        gene_trees[[i]]$tip.label,
        gene_trees[[j]]$tip.label
      )
      
      if (length(common_tips) < min_tips) next
      
      ti <- drop.tip(
        gene_trees[[i]],
        setdiff(gene_trees[[i]]$tip.label, common_tips)
      )
      tj <- drop.tip(
        gene_trees[[j]],
        setdiff(gene_trees[[j]]$tip.label, common_tips)
      )
      
      rf_vec <- c(
        rf_vec,
        RF.dist(
          ti, tj,
          normalize = normalize,
          rooted = rooted
        )
      )
    }
  }
  
  rf_mat <- matrix(rf_vec, ncol = 1)
  colnames(rf_mat) <- "RF"
  
  return(rf_mat)
}

# -------------------------------------------------------------------------

## 设置树的路径
all_tree <- list.files("./tree/orthofinder_genes/",pattern = "\\.tre$",full.names = T)


## 4.1 按照分布密度来计算------------------------------------------------------
## 区间模式选择
## "paired"     : 隔一个取一段 (1-2, 3-4, 5-6, ...)
## "consecutive": 取所有连续段 (1-2, 2-3, 3-4, 4-5, ...)
interval_mode <- "paired"
# 设置分位点（成对出现，每对定义一个区间）
# orthofinder_paired
# quantile_points <- c(0.40, 0.45, 0.5, 0.55, 0.605, 0.65,0.71,0.76,0.82,0.86)
quantile_points <- c(0.40, 0.45, 0.5, 0.55, 0.605, 0.65,0.71,0.76,0.82,0.86)
n_intervals=5
# 设置每个区间的颜色
interval_colors <- c("grey","#4DAF4A", "#377EB8","orange", "#E41A1C")

# 设置每个区间的标签（可选）
interval_labels <- c("Low","Low_m","Mid","Mid_h","High")

## 主程序 
merged <- merged_orthofinder
factor <- "Proportion_parsimony_informative"

# ================================================

# 计算密度
d <- density(merged[[factor]], adjust = 3)
df <- data.frame(x = d$x, y = d$y)

# 计算实际分位点值
x <- merged[[factor]]
qs <- round(
  min(x, na.rm = TRUE) + 
    quantile_points * (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
  2
)


interval_data <- list()
genes_list <- list()

## 按照数值比例区间填充“interval_data”和“genes_list”
for (i in 1:n_intervals) {
  # 根据模式确定起止索引
  if (interval_mode == "paired") {
    idx_start <- 2 * i - 1
    idx_end <- 2 * i
  } else {  # consecutive
    idx_start <- i
    idx_end <- i + 1
  }
  
  # 提取区间数据
  interval_data[[i]] <- df[df$x >= qs[idx_start] & df$x <= qs[idx_end], ]
  
  # 提取基因名
  genes_list[[i]] <- merged[["Alignment_name"]][
    merged[[factor]] >= qs[idx_start] & 
      merged[[factor]] <= qs[idx_end]
  ]
  
  # 设置名称
  if (length(interval_labels) >= i) {
    names(genes_list)[i] <- interval_labels[i]
  } else {
    names(genes_list)[i] <- paste0("Interval_", i)
  }
}



## 绘图：使用循环避免重复代码 ---
# 初始化基础图层
p1 <- ggplot(df, aes(x, y)) +
  geom_area(fill = "#7393B3", alpha = 0.75) +
  geom_line(linewidth = 0.6) +
  theme_bw() +
  labs(x = factor, y = "Density")

# 循环添加区间和标签
for (i in 1:n_intervals) {
  idx_start <- 2 * i - 1
  idx_end <- 2 * i
  
  # 添加区间 ribbon
  p1 <- p1 + 
    geom_ribbon(
      data = interval_data[[i]],
      aes(ymin = 0, ymax = y),
      fill = interval_colors[i],
      alpha = 0.6
    )

  
  
  # 添加起点和终点标签
  p1 <- p1 + 
    annotate("text", x = qs[idx_start], y = 0, 
             label = qs[idx_start],
             vjust = 1.5, size = 3) +
    annotate("text", x = qs[idx_end], y = 0, 
             label = qs[idx_end],
             vjust = 1.5, size = 3)
}


# 显示图形
print(p1)

df_pvs <- do.call(rbind, lapply(1:n_intervals, function(i) {
  data.frame(
    Proportion_variable_sites = merged[[factor]][
      merged[["Alignment_name"]] %in% genes_list[[i]]
    ],
    Group = names(genes_list)[i]
  )
}))

df_pvs$Group <- factor(df_pvs$Group, levels = names(genes_list))

p2 <- ggplot(df_pvs, aes(x = Group, y = Proportion_variable_sites, fill = Group)) +
  geom_violin(
    width = 1,
    alpha = 0.4,
    trim = TRUE
  ) +
  geom_boxplot(
    fill="white",
    width = 0.2,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = interval_colors
  ) +
  theme_classic(base_size = 14) +
  ylab(factor) +
  xlab("")+
  theme(legend.position = "none")

## 4.2 不使用数值区间，而是按照一定的基因数量直接等距划分基因集------------------
interval_colors <- c("grey","#4DAF4A", "#377EB8","orange", "#E41A1C","black")

# 设置每个区间的标签（可选）
interval_labels <- c("Low","Low_m","Mid","Mid_h","High")
# 设置梯度
quantile_points <- list(
  Low  = c(0.1, 0.15),
  Low_m  = c(0.25, 0.3),
  Mid = c(0.4, 0.45),
  Mid_h = c(0.55,0.6),
  High = c(0.7,0.75)
)

# 按照分析的因子对数据框进行重排序
merged <- merged[order(merged[[factor]]),]
n_total <- nrow(merged)

# 初始化储存对象
extracted_data <- data.frame()
genes_list <- list()

#循环提取所需数据
df_pvs <- data.frame()
for(group_name in names(quantile_points)){
  start <- quantile_points[[group_name]][1]
  end <- quantile_points[[group_name]][2]
  start_index <- floor(n_total*start)+1
  end_index <- floor(n_total*end)
  tmp_df <- merged[start_index:end_index,c("Alignment_name",factor)]
  tmp_df$Group <- group_name
  df_pvs <- rbind(df_pvs,tmp_df)
  genes_list[[group_name]] <- merged$Alignment_name[start_index:end_index]
}


# 为原始数据打上分组标签
merged$Group <- "Others"
for (group_name in names(quantile_points)) {
  start_pct <- quantile_points[[group_name]][1]
  end_pct   <- quantile_points[[group_name]][2]
  
  start_idx <- floor(nrow(merged) * start_pct) + 1
  end_idx   <- floor(nrow(merged) * end_pct)
  
  merged$Group[start_idx:end_idx] <- group_name
}
# 设置因子水平，确保图例顺序为 Low -> Mid -> High
merged$Group <- factor(merged$Group, levels = c(names(quantile_points), "Others"))



p1 <- ggplot(merged, aes(x = Rank, y = .data[[factor]], fill = Group)) +
  geom_col(width = 0.8) +  
  scale_fill_manual(
    values = interval_colors,
    # 直接使用 names(intervals) 作为显示的断点，排除 Others
    breaks = names(quantile_points) 
  ) +
  geom_vline(
    xintercept = sapply(quantile_points, function(x) floor(nrow(merged) * x[2])) + 0.5,
    linetype = "dashed", 
    color = "black", 
    alpha = 0.3
  ) +
  labs(
    title = paste("Data Partitioning based on", factor),
    x = "Genes (Ranked)",
    y = factor
  ) +
  theme_minimal() +
  theme(legend.position = "top")
print(p1)


df_pvs$Group <- factor(df_pvs$Group, levels = c(names(quantile_points)))
p2 <- ggplot(df_pvs, aes(x = Group, y = Proportion_parsimony_informative, fill = Group)) +
  geom_violin(
    width = 1,
    alpha = 0.4,
    trim = TRUE
  ) +
  geom_boxplot(
    fill="white",
    width = 0.2,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = interval_colors
  ) +
  theme_classic(base_size = 14) +
  ylab(factor) +
  xlab("")+
  theme(legend.position = "none")

print(p2)




## 4.3开始读取所有的基因树，并计算RF距离---------------------------

### -----------------------------------------------------------

# 接4.2 正式计算
rf_list <- list()
rf_list <- lapply(1:n_intervals, function(i) {
  pairwise_rf_onecol(
    genes = genes_list[[i]],
    pattern_suffix = "_",
    all_tree,
    min_tips = 5
  )
})

# 动态构建 RF 数据框
df_rf <- do.call(rbind, lapply(1:n_intervals, function(i) {
  data.frame(
    RF = as.numeric(rf_list[[i]]),
    Group = names(genes_list)[i]
  )
}))

# 接4.3
rf_list <- lapply(names(quantile_points), function(group_name) {
  
  message(paste("正在计算组别:", group_name, "的 RF 距离..."))
  
  # 从之前生成的 genes_list 中提取对应组的基因名
  current_genes <- genes_list[[group_name]]
  
  # 调用你的函数
  # 注意：确保 pairwise_rf_onecol 能处理这些基因名
  res <- pairwise_rf_onecol(
    genes = current_genes,
    pattern_suffix = "_",
    all_tree = all_tree, # 确保变量名一致
    min_tips = 5
  )
  
  return(res)
})

names(rf_list) <- names(quantile_points)

# 动态构建 RF 数据框
df_rf <- do.call(rbind, lapply(names(quantile_points), function(group_name) {
  data.frame(
    RF = as.numeric(rf_list[[group_name]]),
    Group = group_name
  )
}))



   
# 设置分组因子顺序
df_rf$Group <- factor(df_rf$Group, levels = names(genes_list))


# 设置分组因子顺序
df_rf$Group <- factor(df_rf$Group, levels = names(genes_list))

gene_n <- sapply(genes_list, length)
group_labels <- paste0(names(gene_n), "\n(", gene_n, ")")
names(group_labels) <- names(gene_n)

p3 <- ggplot(df_rf, aes(x = Group, y = RF, fill = Group)) +
  geom_violin(
    width = 1,
    alpha = 0.4,
    trim = TRUE
  ) +
  geom_boxplot(
    fill="white",
    width = 0.2,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = interval_colors
  ) +
  scale_x_discrete(labels = group_labels) +
  theme_classic(base_size = 14) +
  ylab("RF") +
  xlab("")+
  theme(legend.position = "none")




pdf("orthofinder_parsimony_rf_equidistant.pdf",width = 12,height = 10)
plot_grid(
  p1,
  plot_grid(
    p2, p3, 
    ncol = 2, 
    align = "h",  # 关键：水平对齐面板
    axis = "tb",  # 关键：对齐顶部和底部坐标轴
    rel_widths = c(1, 1)
  ),
  ncol = 1,
  rel_widths = c(1, 1)
)
dev.off()

rm(col,d,df,df_long,df_pvs,df_rf,final,genes_list,high,interval_data,low,low_m,mid,mid_h,rf_list,rf_high,rf_mid,rf_low,res,p_mat,cor_mat)


# 5.寻找orthofinder与chloroplast相似的基因 (完成此步后，需建树)----------------------------------------------------------------

## 5.1 建立函数 -----------------------------------------------------------------
## 寻找与质体最匹配的一个基因---自定义函数 
match_property_once <- function(chloro, nuc, prop, pct=0.25){
  
  chloro <- chloro[order(chloro[[prop]]), ]
  nuc    <- nuc[order(nuc[[prop]]), ]
  
  selected <- c()
  
  result <- data.frame(
    chloro_gene   = character(),
    chloro_value  = numeric(),
    matched_gene  = character(),
    matched_value = numeric(),
    diff          = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(i in seq_len(nrow(chloro))){
    
    target_val <- chloro[[prop]][i]
    lower <- target_val  - abs(target_val*pct)
    upper <- target_val + abs(target_val*pct)
    
    cand_idx <- which(
      nuc[[prop]] >= lower &
        nuc[[prop]] <= upper &
        !(nuc$Alignment_name %in% selected)
    )
    
    if(length(cand_idx)==0){
      matched_gene  <- NA
      matched_value <- NA
      diff <- NA
      
    } else {
      
      diffs <- abs(nuc[[prop]][cand_idx] - target_val)
      picked_idx <- cand_idx[which.min(diffs)]
      
      matched_gene  <- nuc$Alignment_name[picked_idx]
      matched_value <- nuc[[prop]][picked_idx]
      diff <- abs(target_val - matched_value)
      selected <- c(selected, matched_gene)
    }
    
    result <- rbind(
      result,
      data.frame(
        chloro_gene   = chloro$Alignment_name[i],
        chloro_value  = target_val,
        matched_gene  = matched_gene,
        matched_value = matched_value,
        diff          = diff
      )
    )
  }
  
  return(result)
}



## 寻找与质体最相近的2个基因，用于后续随机抽样---自定义函数
match_property_once_random <- function(chloro, nuc, prop, pct=0.25){
  
  chloro <- chloro[order(chloro[[prop]]), ]
  nuc    <- nuc[order(nuc[[prop]]), ]
  
  result <- data.frame(
    chloro_gene   = character(),
    chloro_value  = numeric(),
    matched_gene  = character(),
    matched_value = numeric(),
    diff          = numeric(),
    all_nuc_gene  = character(),
    n_matched     = integer(),
    stringsAsFactors = FALSE
  )
  
  # === 第一轮：为每个叶绿体基因匹配最相似的1个核基因 ===
  selected <- c()  # 已被选中的核基因索引
  first_match <- list()  # 存储第一轮匹配结果
  
  for(i in seq_len(nrow(chloro))){
    target_val <- chloro[[prop]][i]
    lower <- target_val  - abs(target_val*pct)
    upper <- target_val + abs(target_val*pct)
    
    # 找到所有候选（未被选中）
    cand_idx <- which(
      nuc[[prop]] >= lower &
        nuc[[prop]] <= upper &
        !(seq_len(nrow(nuc)) %in% selected)
    )
    
    if(length(cand_idx) == 0){
      first_match[[i]] <- list(idx = NULL, gene = NA, value = NA, diff = NA)
    } else {
      # 找最相似的1个
      diffs <- abs(target_val - nuc[[prop]][cand_idx])
      best_idx <- cand_idx[which.min(diffs)]
      
      first_match[[i]] <- list(
        idx   = best_idx,
        gene  = nuc$Alignment_name[best_idx],
        value = nuc[[prop]][best_idx],
        diff  = abs(target_val - nuc[[prop]][best_idx])
      )
      
      # 标记为已选中
      selected <- c(selected, best_idx)
    }
  }
  
  # === 第二轮：为每个叶绿体基因匹配第2个最相似的核基因 ===
  second_match <- list()
  
  for(i in seq_len(nrow(chloro))){
    target_val <- chloro[[prop]][i]
    lower <- target_val  - abs(target_val*pct)
    upper <- target_val + abs(target_val*pct)
    
    # 找到所有候选（未被选中，且不是自己第一轮选的）
    cand_idx <- which(
      nuc[[prop]] >= lower &
        nuc[[prop]] <= upper &
        !(seq_len(nrow(nuc)) %in% selected)
    )
    
    if(length(cand_idx) == 0){
      second_match[[i]] <- list(idx = NULL, gene = NA)
    } else {
      # 找最相似的1个
      diffs <- abs(target_val - nuc[[prop]][cand_idx])
      best_idx <- cand_idx[which.min(diffs)]
      
      second_match[[i]] <- list(
        idx  = best_idx,
        gene = nuc$Alignment_name[best_idx]
      )
      
      # 标记为已选中
      selected <- c(selected, best_idx)
    }
  }
  
  # === 整合结果 ===
  for(i in seq_len(nrow(chloro))){
    first  <- first_match[[i]]
    second <- second_match[[i]]
    
    # 收集所有匹配的基因名
    all_genes <- c()
    if(!is.na(first$gene)) all_genes <- c(all_genes, first$gene)
    if(!is.na(second$gene)) all_genes <- c(all_genes, second$gene)
    
    result <- rbind(
      result,
      data.frame(
        chloro_gene   = chloro$Alignment_name[i],
        chloro_value  = chloro[[prop]][i],
        matched_gene  = first$gene,
        matched_value = first$value,
        diff          = first$diff,
        all_nuc_gene  = if(length(all_genes) > 0) paste(all_genes, collapse = "; ") else NA,
        n_matched     = length(all_genes)
      )
    )
  }
  
  return(result)
}


## 从上一个函数的结果出发，随机抽样20次---自定义函数
random_sample_from_candidates <- function(candidates_df, n_iterations=20){
  
  all_results <- list()
  
  for(iter in 1:n_iterations){
    
    iter_result <- data.frame(
      chloro_gene   = character(),
      matched_gene  = character(),
      iteration     = integer(),
      stringsAsFactors = FALSE
    )
    
    for(i in seq_len(nrow(candidates_df))){
      
      chloro_gene <- candidates_df$chloro_gene[i]
      all_nuc <- candidates_df$all_nuc_gene[i]
      
      if(is.na(all_nuc)){
        # 没有候选，记录NA
        matched_gene <- NA
      } else {
        # 分割候选基因
        nuc_genes <- strsplit(all_nuc, "; ")[[1]]
        
        # 随机抽取1个
        matched_gene <- sample(nuc_genes, 1)
      }
      
      iter_result <- rbind(
        iter_result,
        data.frame(
          chloro_gene  = chloro_gene,
          matched_gene = matched_gene,
          iteration    = iter,
          stringsAsFactors = FALSE
        )
      )
    }
    
    all_results[[iter]] <- iter_result
  }
  
  # 合并所有迭代结果
  final_result <- do.call(rbind, all_results)
  
  return(final_result)
}

## 绘制20次重复的密度分布图---自定义函数
plot_random_iterations <- function(random_results, nuc_data, chloro_data, all_nuc_data, prop = "Mean_evolutionary_rate"){

  # 创建基础图层
  p <- ggplot() +
    # 背景：所有核基因
    geom_density(data = all_nuc_data, aes(x = .data[[prop]]),
                 adjust = 3, fill = "yellow", alpha = 0.6) +
    # 叶绿体基因
    geom_density(data = chloro_data, aes(x = .data[[prop]]),
                 adjust = 3, fill = "#a64d6d", alpha = 0.6)
  
  # 为每次迭代添加密度曲线
  for(iter in unique(random_results$iteration)){
    
    # 提取该迭代的匹配基因
    iter_genes <- random_results %>%
      filter(iteration == iter, !is.na(matched_gene)) %>%
      pull(matched_gene)
    
    # 从核基因数据中提取对应的属性值
    iter_data <- nuc_data %>%
      filter(Alignment_name %in% iter_genes)
    
    # 添加密度曲线
    if(nrow(iter_data) > 0){
      p <- p + geom_density(data = iter_data, aes(x = .data[[prop]]),
                            adjust = 3, fill = "#2073a5", alpha = 0.1)
    }
  }
  
  # 添加标签和主题
  p <- p +
    labs(x = prop, y = "Density", 
         title = paste0("Distribution of ", prop, " across 20 random iterations")) +
    theme_bw()
  
  return(p)
}


## 对叶绿体进行抽样稀疏,按照每3个抽1个的规则进行抽样---自定义函数
sample_every_n_genes <- function(chloro_data, prop = "Mean_evolutionary_rate", interval = 3){
  
  library(dplyr)
  
  # 按属性值排序
  chloro_sorted <- chloro_data %>%
    arrange(.data[[prop]])
  
  n_total <- nrow(chloro_sorted)
  
  # 每interval个抽1个
  selected_indices <- seq(1, n_total, by = interval)
  
  # 如果想从中间位置开始抽（比如每3个中抽第2个），可以调整起始位置
  # selected_indices <- seq(2, n_total, by = interval)  # 从第2个开始
  
  result <- chloro_sorted[selected_indices, ]
  
  # 打印抽样信息
  cat("Sampling summary:\n")
  cat(sprintf("  Total genes in dataset: %d\n", n_total))
  cat(sprintf("  Sampling interval: every %d genes\n", interval))
  cat(sprintf("  Sampled genes: %d\n", nrow(result)))
  cat(sprintf("  GC_content range: %.4f - %.4f\n", 
              min(result[[prop]]), max(result[[prop]])))
  
  return(result)
}


# 
# # 使用示例
# sampled_chloroplast <- sample_every_n_genes(
#   merged_chloroplast, 
#   prop = "Mean_evolutionary_rate", 
#   interval = 3)



## 5.2 所有函数加载完毕，找出Chloroplast-like核基因并绘图------------------------
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"GC_content", pct=0.35)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"Proportion_variable_sites", pct=0.45)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"Proportion_parsimony_informative", pct=0.75)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"RCV", pct=2.4)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"PC1", pct=0.55)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"Mean_evolutionary_rate_without0", pct=2)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"Mean_internal_branch", pct=8)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,"Mean_teminal_branch", pct=3.5)



## 开始匹配，配置factor,并注意调控第一个match_property_once中的pct值

factor <- "PC1"

merged_chloroplast[[factor]] <-
  as.numeric(as.character(merged_chloroplast[[factor]]))

merged_orthofinder[[factor]] <-
  as.numeric(as.character(merged_orthofinder[[factor]]))

# 调控pct值，使得至少能够匹配出与叶绿体基因树量相同的核基因
# 在此前提下，pct越小越好，越大能匹配出的核基因数量越多，但同时也越不精准
# 同时也是由于pct值需要调试，目前难以直接批量化循环
result_gene_GC <- match_property_once(merged_chloroplast, merged_orthofinder,factor, pct=0.55)

result_gene_GC <- match_property_once_random(merged_chloroplast, merged_orthofinder,factor, pct=0.55)

random_results <- random_sample_from_candidates(result_gene_GC,n_iterations=20)



## 使用示例
p <- plot_random_iterations(
  random_results = random_results,
  nuc_data = merged_orthofinder,
  chloro_data = merged_chloroplast,
  all_nuc_data = merged_orthofinder,
  prop = factor
)

print(p)


## 随机抽取20个核基因作为重复

## 设置参数
n_iterations <- 20
sample_size <- 75

## 20次随机抽样，每次画一条黄色密度曲线
set.seed(123)
sampling_results <- list()
for(i in 1:n_iterations){
  sampled <- merged_orthofinder[sample(nrow(merged_orthofinder), sample_size), ]
  # 添加迭代标记
  sampled_with_iter <- sampled %>%
    select(Alignment_name, factor) %>%
    mutate(Iteration = i)
  
  # 存储到列表
  sampling_results[[i]] <- sampled_with_iter
  
  p <- p + geom_density(data = sampled, aes(x = .data[[factor]]),
                        adjust = 3, fill = "yellow", alpha = 0.1)
}

## 添加标签和主题
p <- p +
  labs(x = factor, y = "Density",
       title = paste0("Distribution of ",factor," ", n_iterations, " random iterations (n=", sample_size, ")")) +
  theme_bw()+
  geom_point(aes(x = -Inf, y = -Inf, color = "Random nuclear genes"), size = 0) +
  geom_point(aes(x = -Inf, y = -Inf, color = "Chloroplast like gene"), size = 0) +
  geom_point(aes(x = -Inf, y = -Inf, color = "Chloroplast genes"), size = 0) +
  scale_color_manual(
    name = "Dataset",
    values = c(
      "Random nuclear genes" = "yellow",
      "Chloroplast like gene" = "#2073a5",
      "Chloroplast genes" = "#a64d6d"
    ),
    labels = c(
      "Random nuclear genes" = "Random nuclear genes",
      "Chloroplast like gene" = "Chloroplast like gene",
      "Chloroplast genes" = "Chloroplast genes"
    )
  )+
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 0.6))) +
  theme(legend.position = c(0.85, 0.85),  # 右上角位置 (x, y)，范围0-1
        legend.background = element_rect(fill = "white", color = "black"),  # 图例背景
        legend.title = element_text(face = "bold"))

print(p)

ggsave(paste0("./distribution_factor_plot_20/",factor,"_random_iterations.pdf"),p, width = 10, height = 6)


## 5.3 将最终挑选出来的20次重复采样的叶绿体相似核基因读出-----------------------
## 合并所有迭代结果为一个数据框
all_sampling_results <- bind_rows(sampling_results)

## 重新排列列顺序
all_sampling_results <- all_sampling_results %>%
  select(Alignment_name, Iteration, factor)

## 转换为宽数据框
all_sampling_results_width <- all_sampling_results %>%
  group_by(Iteration) %>%
  mutate(row_id = row_number()) %>%  # 为每次迭代内的基因编号
  ungroup() %>%
  pivot_wider(
    id_cols = row_id,
    names_from = Iteration,
    values_from = Alignment_name,
    names_prefix = "Iteration_"
  ) %>%
  select(-row_id)

random_results_width <- random_results %>% 
  group_by(iteration) %>%
  mutate(row_id = row_number()) %>%  # 为每次迭代内的基因编号
  ungroup() %>%
  pivot_wider(
    id_cols = row_id,
    names_from = iteration,
    values_from = matched_gene,
    names_prefix = "Iteration_"
  ) %>%
  select(-row_id)

write.csv(random_results,paste0("./distribution_factor_plot_20/",factor,"_Chloroplast_like_nuclear_gene.csv"),row.names = F,quote = F)
write.csv(random_results_width,paste0("./distribution_factor_plot_20/",factor,"_Chloroplast_like_nuclear_gene_genelist.csv"),row.names = F,quote = F)



# 6.比较不同的因素中的like-chloroplast核基因树与质体树的RF距离 --------------------------------------------------------

## 自定义函数计算叶绿体树与质体基因间的差异
get_rf_distanch <- function(tree_list,chloroplast,name){
  rf_list <- numeric(length(tree_list))
  for (i in seq_along(tree_list)){
    tt <- read.tree(tree_list[i])
    rf_list[i] <- RF.dist(tt, chloroplast, check.labels = TRUE, normalize = TRUE)
  }
  df_chloroplast <- data.frame(rf = rf_list)
  df_chloroplast$type <- name
  return(df_chloroplast)
}


## 读取质体树
chloroplast <- read.tree("tree/rosa_chloroplast_partition_rt.tre")

## 随机基因树
tree_list <- list.files("./tree/Random_nuclear/",pattern = "\\.tre$",full.names = T)
df_random <- get_rf_distanch(tree_list,chloroplast,"Random_nuclear")

## GC含量
tree_list <- list.files("./tree/GC_content/",pattern = "\\.tre$",full.names = T)
df_GC <- get_rf_distanch(tree_list,chloroplast,"GC_content")

## Mean_evolutionary_rate
tree_list <- list.files("./tree/Mean_evolutionary_rate/",pattern = "\\.tre$",full.names = T)
df_rate <- get_rf_distanch(tree_list,chloroplast,"Mean_evolutionary_rate")

## PC1
tree_list <- list.files("./tree/PC1/",pattern = "\\.tre$",full.names = T)
df_PC1 <- get_rf_distanch(tree_list,chloroplast,"PC1")

## Proportion_parsimony_informative
tree_list <- list.files("./tree/Proportion_parsimony_informative/",pattern = "\\.tre$",full.names = T)
df_parsimony <- get_rf_distanch(tree_list,chloroplast,"Proportion_parsimony_informative")

## Proportion_variable_sites
tree_list <- list.files("./tree/Proportion_variable_sites/",pattern = "\\.tre$",full.names = T)
df_var <- get_rf_distanch(tree_list,chloroplast,"Proportion_variable_sites")

## RCV
tree_list <- list.files("./tree/RCV/",pattern = "\\.tre$",full.names = T)
df_RCV <- get_rf_distanch(tree_list,chloroplast,"RCV")

## Mean internal branch
tree_list <- list.files("./tree/Mean_internal_branch/",pattern = "\\.tre$",full.names = T)
df_internal <- get_rf_distanch(tree_list,chloroplast,"Mean_internal_branch")

## Mean teminal branch
tree_list <- list.files("./tree/Mean_teminal_branch/",pattern = "\\.tre$",full.names = T)
df_teminal <- get_rf_distanch(tree_list,chloroplast,"Mean_teminal_branch")



## 画图展示差异
df <- bind_rows(df_random,df_GC,df_rate,df_PC1,df_parsimony,df_var,df_RCV,df_internal,df_teminal)
df$type <- factor(
  df$type,
  levels = c(
    "Random_nuclear",
    "GC_content",
    "Mean_evolutionary_rate",
    "PC1",
    "Proportion_parsimony_informative",
    "Proportion_variable_sites",
    "RCV",
    "Mean_internal_branch",
    "Mean_teminal_branch"
  )
)


## 提取 baseline

med=median(df_random$rf)

p <- ggplot(df, aes(x = type, y = rf, fill = type)) +
  geom_boxplot() +
  geom_hline(yintercept=med, color="red", linetype="dashed")+
  scale_fill_manual(
    values = c(
      "Random_nuclear" = "white",
      "GC_content" = "grey70",
      "Mean_evolutionary_rate" = "grey70",
      "PC1" = "grey70",
      "Proportion_parsimony_informative" = "grey70",
      "Proportion_variable_sites" = "grey70",
      "RCV" = "grey70",
      "Mean_internal_branch" = "grey70",
      "Mean_teminal_branch" = "grey70"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none"
  )

ggsave("Boxplot_all_diff.pdf",p, width = 20, height = 10)


