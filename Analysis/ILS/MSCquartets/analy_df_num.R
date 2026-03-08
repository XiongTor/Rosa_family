#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2026-02-06
# Description: Analyse the results from the 'results_classify.sh' file.

# ==== Main Script Start ====
library(ggplot2)
library(ggrepel) # 用于自动标注离群基因的标签
library(dplyr)
library(tidyverse)
library(UpSetR)

# 1. 读入数据
df_T3 <- read.csv("./T3/all_df_num.csv")
df_T1 <- read.csv("./T1/all_df_num.csv")

data <- c("df_T1","df_T3")
alpha_list <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001)
beta_list <- c(0.01, 0.05)

all_outlier_red <- data.frame()
all_outlier_yellow <- data.frame()


for(name in data){
  for (alpha in alpha_list) {
    for (beta in beta_list) {
      df <- get(name)
      df <- df %>%
        mutate(
          # 权重系数：利用 log10 缩放四分体总数
          # 为了防止 log10(1) = 0 导致结果消失，可以统一加 1
          weight = log10(all_quartet + 1),
          
          # 计算加权指数：IH 概率 * 权重
          weighted_red = red * weight,
          weighted_yellow=yellow * weight
        )
            
      # 按照不同的alpha和beta值取样
      df1 <- df %>% 
        filter(df$alpha_value== alpha & df$beta_value==beta)
      
      # # 四分法计算离群值
      # # # 2. 计算 Red 与Yellow的阈值，用于识别离群基因
      # # # Red，计算 IQR 阈值
      # Q1 <- quantile(df1$weighted_red, 0.25)
      # Q3 <- quantile(df1$weighted_red, 0.75)
      # IQR_val <- Q3 - Q1
      # red_threshold <- Q3 + 1.5 * IQR_val
      # 
      # # Yellow,
      # Q1 <- quantile(df1$weighted_yellow, 0.25)
      # Q3 <- quantile(df1$weighted_yellow, 0.75)
      # IQR_val <- Q3 - Q1
      # yellow_threshold <- Q3 + 1.5 * IQR_val
      
      # 标准差计算离群值
      red_threshold <- mean(df1$weighted_red) + 2 * sd(df1$weighted_red)
      yellow_threshold <- mean(df1$weighted_yellow) + 2 * sd(df1$weighted_yellow)
      # 3. 标记哪些基因是显著的渗入候选
      df1$label <- ifelse(df1$weighted_red > red_threshold | df1$weighted_yellow > yellow_threshold, as.character(df1$gene), "")
      

      
      # 4. 开始绘图
      p <- ggplot(df1, aes(x = weighted_yellow, y = weighted_red)) +
        # 绘制所有基因的散点
        geom_point(aes(color = red), alpha = 0.6) +
        # 映射颜色梯度（学术风常用蓝色到红色）
        scale_color_gradient(low = "#3C5488", high = "#E64B35") +
        # 绘制红色阈值虚线
        geom_hline(yintercept = red_threshold, linetype = "dashed", color = "red") +
        geom_vline(xintercept = yellow_threshold, linetype = "dashed", color = "red") +
        # 标注 Red 占比极高的离群基因 ID
        geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
        theme_bw() +
        labs(
          title = "All genes contribution for IH/ILS",
          subtitle = paste0("alpha=",alpha," ","beta=",beta),
          x = "Yellow (ILS)",
          y = "Red (Hybridization)",
          color = "Red Score"
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          panel.grid.minor = element_blank()
        )
      
      ggsave(paste0("All_gene_contribution","_alpha",alpha,"_beta",beta,"_",name,".pdf"),p,width = 10,height = 8)
      
      # 提取具体的基因名称 -------------------------------------------------------
      red_gene<- df1[df1$weighted_red>red_threshold,]
      write.csv(red_gene,paste0("Red_gene","_alpha",alpha,"_beta",beta,"_",name,".csv"))
      
      yellow_gene <-  df1[df1$weighted_yellow>yellow_threshold,]
      write.csv(yellow_gene,paste0("Yellow_gene","_alpha",alpha,"_beta",beta,"_",name,".csv"))
      
      # 全部结果合并到一个数据集中
      if(nrow(red_gene) > 0) {
        all_outlier_red <- rbind(all_outlier_red, red_gene)
      }
      
      if(nrow(yellow_gene) > 0) {
        all_outlier_yellow <- rbind(all_outlier_yellow, yellow_gene)
      }
    }
  }
}

all_outlier_red_unique <- all_outlier_red %>% distinct()
all_outlier_yellow_unique <- all_outlier_yellow %>% distinct()


write.csv(all_outlier_red_unique,"all_outlier_red.csv",quote = F,row.names = F)
write.csv(all_outlier_red_unique,"all_outlier_yellow.csv",quote = F,row.names = F)




# # 四分法计算离群值
# # # 2. 计算 Red 与Yellow的阈值，用于识别离群基因
# # # Red，计算 IQR 阈值
# Q1 <- quantile(df$weighted_red, 0.25)
# Q3 <- quantile(df$weighted_red, 0.75)
# IQR_val <- Q3 - Q1
# red_threshold <- Q3 + 1.5 * IQR_val
# 
# # Yellow,
# Q1 <- quantile(df$weighted_yellow, 0.25)
# Q3 <- quantile(df$weighted_yellow, 0.75)
# IQR_val <- Q3 - Q1
# yellow_threshold <- Q3 + 1.5 * IQR_val

# -----------------------------------------------------------------------------
# 统计所有alpha和beta组合中的基因的并集
# 统计交集
# 展示不同组合下，具体的gene是如何增减的

red_data <- all_outlier_red_unique %>% 
  mutate(comb=paste0("alpha",alpha_value,"beta",beta_value)) %>% 
  select(gene,comb) %>% 
  distinct() 

# 1. 交集：所有条件下共有的基因
gene_list <- red_data %>%
  split(.$comb) %>%
  lapply(function(x) as.character(x$gene))

common_genes <- Reduce(intersect, gene_list)

cat("所有条件下共有的基因数量:", length(common_genes), "\n")

write.csv(data.frame(gene = common_genes), "Common_genes_all_conditions.csv",row.names = FALSE)

# 2. 并集：至少在一个条件下出现的所有基因
all_genes <- Reduce(union, gene_list)

cat("至少出现一次的基因总数:", length(all_genes), "\n")

write.csv(data.frame(gene = all_genes),"Union_genes_all_conditions.csv",row.names = FALSE)

# 3. 热图展示变化
gene_comb_matrix <- red_data %>%
  mutate(present = 1) %>%
  complete(gene, comb, fill = list(present = 0)) %>%
  pivot_wider(names_from = comb, 
              values_from = present, 
              values_fill = 0,
              id_cols = gene)

# 计算每个基因出现的次数，用于排序
gene_order <- gene_comb_matrix %>%
  pivot_longer(-gene, names_to = "comb", values_to = "present") %>%
  group_by(gene) %>%
  summarise(total = sum(present)) %>%
  arrange(desc(total), gene)

# 转换为长格式用于绘图
plot_data <- gene_comb_matrix %>%
  pivot_longer(-gene, names_to = "comb", values_to = "present") %>%
  mutate(gene = factor(gene, levels = gene_order$gene))

# 按照alpha和beta值排序comb
plot_data <- plot_data %>%
  mutate(comb = factor(comb, levels = c("alpha0.01beta0.01","alpha0.01beta0.05","alpha0.001beta0.01","alpha0.001beta0.05","alpha1e-04beta0.01","alpha1e-04beta0.05","alpha1e-05beta0.01","alpha1e-05beta0.05","alpha1e-06beta0.01","alpha1e-06beta0.05")))

# 绘制热图
p <- ggplot(plot_data, aes(x = comb, y = gene, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.05) +
  scale_fill_manual(values = c("0" = "#F0F0F0", "1" = "#E64B35"),
                    labels = c("Absent", "Present"),
                    name = "Gene Status") +
  theme_minimal() +
  labs(title = "Gene Presence/Absence Pattern Across Conditions",
       subtitle = paste0("Total genes: ", length(unique(plot_data$gene)), 
                         " | Conditions: ", length(unique(plot_data$comb))),
       x = "Parameter Combination (alpha & beta)",
       y = "Gene ID") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("Gene_change_alpha_beta.pdf",p,width = 12,height = 10)