# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.05.04
# Usage: This script is used to count the quartet frequencies from gene trees.

# ==== Main Script Start ====

#install.packages("MSCquartets")

library(ape)
library(Quartet)
library(dplyr)
library(ggplot2)

files <- commandArgs(T)

message("Reading species & gene trees...")
gtrees = read.tree(files[1])
simulate_list = lapply(files[2:5], read.tree)
names(simulate_list) = paste0("sim", 1:4)

### 函数：生成 QT_df + comb
get_qt_df <- function(tree){
  QT = quartetTable(tree, progressbar = TRUE)
  df = as.data.frame(QT)
  df$comb <- apply(df, 1, function(x){
    taxa <- names(x)[which(x == 1)]
    paste(taxa, collapse = " ")
  })
  return(df)
}

message("Calculating quartet tables...")
QT_real   <- get_qt_df(gtrees)
QT_sim_df_list <- lapply(simulate_list, get_qt_df)

# 给 QT_sim_df_list 列加后缀
for(nm in names(QT_sim_df_list)){
  colnames(QT_sim_df_list[[nm]])[1:3] <- paste0(colnames(QT_sim_df_list[[nm]])[1:3], "_", nm)
}

# merge all
merged_df <- QT_real
for(nm in names(QT_sim_df_list)){
  sim_df <- QT_sim_df_list[[nm]]
  merged_df <- merged_df %>%
    left_join(sim_df[,c("comb",grep(nm, names(sim_df), value=T))], by="comb")
}

# 计算比例
calc_prop <- function(df, prefix){
  cols <- paste0(c("12|34","13|24","14|23"), "_", prefix)
  df[,paste0(prefix,"_prop")] <- rowSums(df[,cols], na.rm=TRUE) / rowSums(df[,cols], na.rm=TRUE)
  return(df)
}

merged_df$real_prop <- rowSums(merged_df[,paste0(c("12|34","13|24","14|23"),"_real")]) /
                       rowSums(merged_df[,paste0(c("12|34","13|24","14|23"),"_real")])

for(nm in names(QT_sim_df_list)){
  merged_df[,paste0(nm,"_prop")] <- rowSums(merged_df[,paste0(c("12|34","13|24","14|23"),"_",nm)]) /
                                    rowSums(merged_df[,paste0(c("12|34","13|24","14|23"),"_",nm)])
}

# 线性模型 + Spearman + 绘图
pdf("quartet_correlation_plots.pdf", width=8, height=6)
par(mfrow=c(2,2))

lm_results <- list()
for(nm in names(QT_sim_df_list)){
  x <- merged_df[[paste0(nm,"_prop")]]
  y <- merged_df$real_prop
  fit <- lm(y ~ x)
  lm_results[[nm]] <- summary(fit)
  
  # Spearman
  cor_s <- cor(x, y, method="spearman")
  
  # ggplot
  df_plot <- data.frame(x=x, y=y)
  p <- ggplot(df_plot, aes(x=x, y=y)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", col="red") +
    ggtitle(paste0("Real vs ", nm, "\nSpearman rho=", round(cor_s,3))) +
    xlab(paste0(nm," proportion")) +
    ylab("Real proportion") +
    theme_minimal()
  
  print(p)
}

dev.off()

write.csv(merged_df, "quartet_compare_results.csv", row.names=FALSE)
message("✅ Done! Results saved to quartet_compare_results.csv and quartet_correlation_plots.pdf")

file <- commandArgs(TRUE)

library(ape)
library(phangorn)
library(MSCquartets)

### Reading the observed trees
gtrees=read.newick(file[1])

### Reading the simulated trees
simulate1=read.newick(file[2])
simulate2=read.newick(file[3])
simulate3=read.newick(file[4])
simulate4=read.newick(file[5])
# 统计基因树中所有可能的quartet的三种拓扑出现次数
QT=quartetTable(gtrees,progressbar=T)

QT_df <- as.data.frame(QT)

QT_df$comb <- apply(QT_df, 1, function(x){
  taxa <- names(x)[which(x == 1)]   # 找到该行值为1的列名
  paste(taxa, collapse = " ")       # 用空格拼接
})
