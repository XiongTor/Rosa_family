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
file_ags <- list.files("./ags353/")

df_ags <- list()
for (name in file_ags) {
  df_ags[[name]] <- read.table(paste0("./ags353/",name),header = T)
  df_ags[[name]]$Alignment_name <- sapply(
    strsplit(df_ags[[name]]$Alignment_name,"\\_"),`[`,1
  )
}

merged_ags353 <- reduce(df_ags, full_join, by="Alignment_name") %>% 
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

# 2.画图展示所有因子在叶绿体和核基因之间的数值差异 -------------------------------------------

## 先给数据框加标签方便区分来源
merged_ags353$dataset <- "ags353"
merged_chloroplast$dataset <- "chloroplast"

## 合并
df <- bind_rows(merged_ags353, merged_chloroplast)

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
  all_values <- c(merged_ags353[[col_name]], merged_chloroplast[[col_name]])
  x_min <- min(all_values, na.rm = TRUE)
  x_max <- max(all_values, na.rm = TRUE)
  
  p <- ggplot() +
    geom_density(data = merged_ags353, aes_string(x = col_name),
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
    dens_max <- max(density(merged_ags353[[col_name]], na.rm=TRUE)$y,
                    density(merged_chloroplast[[col_name]], na.rm=TRUE)$y)
    p <- p +
      annotate("rect", xmin = x_min, xmax = x_min + (x_max-x_min)*0.08,
               ymin = dens_max*0.85, ymax = dens_max*0.95, fill="#2073a5") +
      annotate("rect", xmin = x_min, xmax = x_min + (x_max-x_min)*0.08,
               ymin = dens_max*0.65, ymax = dens_max*0.75, fill="#a64d6d") +
      annotate("text", x = x_min + (x_max-x_min)*0.1,
               y = dens_max*0.9, label = "AGS353", hjust = 0, size = 5) +
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

# 3.相关性分析，判断哪些因素相关性高 -------------------------------------------------------------

## 根据上一步的结果，挑选出叶绿体和核基因有差异的因子
diff <- c("GC_content",
  "Proportion_variable_sites","Proportion_parsimony_informative",
  "Mean_evolutionary_rate","RCV")

## 计算相关性
df <- merged_ags353[, diff]

res <- rcorr(as.matrix(df))

cor_mat  <- res$r   # 相关性
p_mat    <- res$P   # p value

pdf("corrplot_diff_factor.pdf",width = 10,height=10)
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
df <- merged_ags353 %>% 
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

# 4.针对上述三个因素挑选高中低三个指标的基因，并计算高中低三个区间内基因树彼此之间的RF距离 ---------------------------------------------------------
merged <- merged_chloroplast
factor <- "Proportion_parsimony_informative"
all_tree <- list.files("./tree/Chloroplast/",pattern = "\\.tre$",full.names = T)

d <- density(merged[[factor]], adjust = 3)
df <- data.frame(x=d$x, y=d$y)
## 计算分位点
x <- merged[[factor]]

qs <- round(
  min(x, na.rm = TRUE) +
    c(0.1, 0.25, 0.35, 0.5, 0.6, 0.75) *
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
  2
)


low <- df[df$x >= qs[[1]] & df$x <= qs[[2]] , ] 
mid <- df[df$x >= qs[[3]] & df$x <= qs[[4]] , ]
high <- df[df$x >= qs[[5]] & df$x <= qs[[6]] , ]

## 提取基因名
genes_low <- merged[["Alignment_name"]][
  merged[[factor]] >= qs[[1]] &
    merged[[factor]] <= qs[[2]]
]

genes_mid <- merged[["Alignment_name"]][
  merged[[factor]] >= qs[[3]] &
    merged[[factor]] <= qs[[4]]
]

genes_high <- merged[["Alignment_name"]][
  merged[[factor]] >= qs[[5]] &
    merged[[factor]] <= qs[[6]]
]



p1 <- ggplot(df, aes(x,y)) +
  # 整体密度
  geom_area(fill="#7393B3", alpha=.75) +
  geom_line(linewidth=0.6) +
  # 只在区间内部画 ribbon（自动贴曲线）
  geom_ribbon(
    data = low,
    aes(ymin=0, ymax=y),
    fill = "#4DAF4A",
    alpha = 0.6
  ) +
  annotate("text", x = qs[[1]], y = 0, label = qs[[1]],
           vjust = 1.5, size=3) +
  annotate("text", x = qs[[2]], y = 0, label = qs[[2]],
           vjust = 1.5, size=3)+
  geom_ribbon(
    data = mid,
    aes(ymin=0, ymax=y),
    fill = "#377EB8",
    alpha = 0.6
  ) +
  annotate("text", x = qs[[3]], y = 0, label = qs[[3]],
           vjust = 1.5, size=3) +
  annotate("text", x = qs[[4]], y = 0, label = qs[[4]],
           vjust = 1.5, size=3)+
  geom_ribbon(
    data = high,
    aes(ymin=0, ymax=y),
    fill = "#E41A1C",
    alpha = 0.6
  ) +
  annotate("text", x = qs[[5]], y = 0, label = qs[[5]],
           vjust = 1.5, size=3) +
  annotate("text", x = qs[[6]], y = 0, label = qs[[6]],
           vjust = 1.5, size=3)+
  labs(x="Proportion_parsimony_informative", y="Density") +
  theme_bw()

p1


## 绘制对应的高中低的因子数值的箱线图
df_pvs <- data.frame(
  Proportion_variable_sites = c(
    merged[[factor]][
      merged[["Alignment_name"]] %in% genes_low
    ],
    merged[[factor]][
      merged[["Alignment_name"]] %in% genes_mid
    ],
    merged[[factor]][
      merged[["Alignment_name"]] %in% genes_high
    ]
  ),
  Group = factor(
    c(
      rep("Low",  length(genes_low)),
      rep("Mid",  length(genes_mid)),
      rep("High", length(genes_high))
    ),
    levels = c("Low", "Mid", "High")
  )
)

low_col  <- "#4DAF4A"  # Low
mid_col  <- "#377EB8"  # Mid
high_col <- "#E41A1C"  # High

p2 <- ggplot(df_pvs, aes(x = Group, y = Proportion_variable_sites, fill = Group)) +
  geom_violin(
    width = 0.6,
    alpha = 0.7,
    trim = TRUE
  ) +
  geom_boxplot(
    fill="white",
    width = 0.1,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c(
      Low  = low_col,
      Mid  = mid_col,
      High = high_col
    )
  ) +
  theme_classic(base_size = 14) +
  ylab(factor) +
  xlab("")+
  theme(legend.position = "none")



## 开始读取所有的基因树，并计算RF距离
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

rf_low <- pairwise_rf_onecol(
  genes = genes_low,
  pattern_suffix = ".",
  all_tree,
  min_tips = 5
)



rf_mid <- pairwise_rf_onecol(
  genes = genes_mid,
  pattern_suffix = ".",
  all_tree,
  min_tips = 5
)

rf_high <- pairwise_rf_onecol(
  genes = genes_high,
  pattern_suffix = ".",
  all_tree,
  min_tips = 5
)

df_rf <- data.frame(
  RF = c(
    as.numeric(rf_low),
    as.numeric(rf_mid),
    as.numeric(rf_high)
  ),
  Group = factor(
    c(
      rep("Low",  nrow(rf_low)),
      rep("Mid",  nrow(rf_mid)),
      rep("High", nrow(rf_high))
    ),
    levels = c("Low", "Mid", "High")
  )
)


p3 <- ggplot(df_rf, aes(x = Group, y = RF, fill = Group)) +
  geom_violin(
    width = 0.6,
    alpha = 0.7,
    trim = TRUE
  ) +
  geom_boxplot(
    fill="white",
    width = 0.1,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = c(
      Low  = low_col,
      Mid  = mid_col,
      High = high_col
    )
  ) +
  theme_classic(base_size = 14) +
  ylab(factor) +
  xlab("")+
  theme(legend.position = "none")


pdf("Chloroplast_parsimony.pdf",width = 10,height = 10)
plot_grid(p2,p3,ncol = 1)
dev.off()




# 5.寻找ags353与chloroplast相似的基因 ----------------------------------------------------------------

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

# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"Proportion_variable_sites", pct=0.9)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"RCV", pct=2.8)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"PC1", pct=1)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"Alignment_length", pct=1)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"Mean_evolutionary_rate_without0", pct=1.5)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"Mean_bipartition", pct=0.8)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"Mean_internal_branch", pct=7)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"Mean_teminal_branch", pct=2.8)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"treeness", pct=0.7)
# result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,"saturation", pct=0.8)


## 开始匹配，配置factor,并注意调控第一个match_property_once中的pct值

factor <- "saturation"

merged_chloroplast[[factor]] <-
  as.numeric(as.character(merged_chloroplast[[factor]]))

merged_ags353[[factor]] <-
  as.numeric(as.character(merged_ags353[[factor]]))


result_gene_GC <- match_property_once(merged_chloroplast, merged_ags353,factor, pct=0.8)

result_gene_GC <- match_property_once_random(merged_chloroplast, merged_ags353,factor, pct=0.8)

random_results <- random_sample_from_candidates(result_gene_GC,n_iterations=20)



## 使用示例
p <- plot_random_iterations(
  random_results = random_results,
  nuc_data = merged_ags353,
  chloro_data = merged_chloroplast,
  all_nuc_data = merged_ags353,
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
  sampled <- merged_ags353[sample(nrow(merged_ags353), sample_size), ]
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
       title = paste0("Distribution of ",factor, n_iterations, " random iterations (n=", sample_size, ")")) +
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

ggsave(paste0(factor,"_random_iterations.pdf"),p, width = 10, height = 6)


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

write.csv(random_results,paste0(factor,"_Chloroplast_like_nuclear_gene.csv"),row.names = F,quote = F)
write.csv(random_results_width,paste0(factor,"_Chloroplast_like_nuclear_gene_genelist.csv"),row.names = F,quote = F)



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



## 画图展示差异
df <- bind_rows(df_random,df_GC,df_rate,df_PC1,df_parsimony,df_var,df_RCV)
df$type <- factor(
  df$type,
  levels = c(
    "Random_nuclear",
    "GC_content",
    "Mean_evolutionary_rate",
    "PC1",
    "Proportion_parsimony_informative",
    "Proportion_variable_sites",
    "RCV"
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
      "RCV" = "grey70"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none"
  )

ggsave("Boxplot_all_diff.pdf",p, width = 20, height = 10)


# 7.PCA -------------------------------------------------------------------------------------------------------------

## 添加基因类型标签
merged_ags353$type <- "nuclear"
merged_chloroplast$type <- "chloroplast"

## 合并
df_all <- bind_rows(merged_ags353, merged_chloroplast)

## 查看有哪些属性列（排除 gene_name 和 type）
setdiff(colnames(df_all), c("Alignment_name", "type"))

## 选属性列（你按实际列名修改）
props <- c("GC_content",
          "Proportion_variable_sites","Proportion_parsimony_informative",
          "Mean_evolutionary_rate","RCV")

## PCA 输入矩阵
X <- scale(df_all[, props])

pca <- prcomp(X, center = TRUE, scale. = TRUE)
saveRDS(df_all, file = "df_all.rds")
saveRDS(pca,    file = "pca.rds")


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
merged_ags353 <- merge(
  merged_ags353,
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


p3 <- ggplot(pca_df, aes(PC2, PC3, color = type)) + 
  geom_point(size = 3, alpha=0.8) +
  stat_ellipse(level = 0.95, size=1) +
  theme_bw(base_size = 16)+
  theme(legend.position = c(0.15,0.9),
        legend.background = element_blank())

plot_grid(p1,p2,p3,ncol=3)
## 三维
library(plotly)

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




