# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.11.14
# This script is used to evaluate the contributions of ILS and IH.

# setwd("2025.11.13")
library(ape)
library(phangorn)
library(MSCquartets)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(doParallel)

message("1. read tree.")
sptree=read.tree("rosa_ags353_treeshrink_sp_rt_oneoutg_final.tre")
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"

message("2. load data")
load("quartet_analysis_results.RData")

#筛选出在不同限制条件下的渐渗物种和高ILS物种
Qtest_2 <- as.data.frame(Qtest_2)

message("3. load function")
##########################   定义若干函数   ##############################
# 简化运算QT
# 原函数
get_qt_df <- function(tree, names){
  QT = quartetTable(tree, names, progressbar = TRUE)
  df = as.data.frame(QT)
  df$comb <- apply(df, 1, function(x){
    taxa <- names(x)[which(x == 1)]
    paste(taxa, collapse = " ")
  })
  return(df)
}


#函数，用于生成comb列
clean_data <- function(df){
  df$comb <- apply(df, 1, function(x){
    taxa <- names(x)[which(x == 1)]
    paste(taxa, collapse = " ")
  })
  df_2 <- df[,(ncol(df)-5):ncol(df)]
  df_2$comb <- sapply(strsplit(df_2$comb, "\\s+"), function(x) {
    paste(x[1:4], collapse = " ")
  })
  return(df_2)
}



#函数，用于判断大小
assign_high_median_low <- function(df, cols){
  # 1. 保留原值
  df2 <- df %>% 
    mutate(across(all_of(cols), ~., .names = "{.col}_value"))
  
  # 2. 行内比较大小并分配 high/median/low
  df2 <- df2 %>% 
    rowwise() %>% 
    mutate(
      ranks = list(rank(-as.numeric(c_across(all_of(cols))), ties.method="first")),
      labels = list(c("high","median","low")[ranks]),
      across(all_of(cols), ~ labels[[which(cols == cur_column())]])
    ) %>% 
    ungroup() %>% 
    select(-ranks, -labels)
  
  return(df2)
}

#函数，用于与IH，ILS数据比对提取
#与IH数据框做比对抓取对应数据,注意其中的df_3是基因树三联体集合中与物种树相违背的三联体数据款
get_contribution <- function(df_3,IH_ILS_level){
  df_b2 <- df_3 %>%
    mutate(
      type = str_extract(comb, "12\\|34|13\\|24|14\\|23"),
      comb = str_trim(str_remove(comb, "12\\|34|13\\|24|14\\|23"))
    )
  
  # 2. 与 df_a 的 comb 进行匹配
  df_m <- df_b2 %>%
    inner_join(IH_ILS_level, by = "comb")
  
  # 3. 根据 type 动态抽取列（含 high/median/low，以及 _value 原始数值）
  mat <- apply(df_m, 2, as.character)
  
  df_m$value <- apply(mat, 1, function(row) {
    t <- row["type"]
    if (is.na(t) || !(t %in% colnames(mat))) return(NA_character_)
    row[t]
  }) 
  
  df_m <- df_m %>% 
    select(c("comb","type","value","12|34_value","13|24_value","14|23_value"))
  
  return(df_m)
}



################   筛选不同等级的IH/ILS的三联体集合   ########################
message("4. count triplet dataframe")
#all
all <- clean_data(Qtest_2)
all<- assign_high_median_low(all,
                             cols = c("12|34","13|24","14|23"))

# 渐渗：
# 0.01
IH_0.01 <- Qtest_2 %>% filter(Qtest_2$p_T3<0.01)
IH_0.01 <- clean_data(IH_0.01)
IH_0.01<- assign_high_median_low(IH_0.01,
                                 cols = c("12|34","13|24","14|23"))


IH_0.001 <- Qtest_2 %>% filter(Qtest_2$p_T3<0.001)
IH_0.001 <- clean_data(IH_0.001)
IH_0.001<- assign_high_median_low(IH_0.001,
                                  cols = c("12|34","13|24","14|23"))

IH_0.0001 <- Qtest_2 %>% filter(Qtest_2$p_T3<0.0001)
IH_0.0001 <- clean_data(IH_0.0001)
IH_0.0001<- assign_high_median_low(IH_0.0001,
                                   cols = c("12|34","13|24","14|23"))

IH_0.00001 <- Qtest_2 %>% filter(Qtest_2$p_T3<0.00001)
IH_0.00001 <- clean_data(IH_0.00001)
IH_0.00001<- assign_high_median_low(IH_0.00001,
                                    cols = c("12|34","13|24","14|23"))


#ILS
ILS_0.05 <- Qtest_2 %>% filter(Qtest_2$p_star>0.05)
ILS_0.05 <- clean_data(ILS_0.05)
ILS_0.05<- assign_high_median_low(ILS_0.05,
                                  cols = c("12|34","13|24","14|23"))

ILS_0.01 <- Qtest_2 %>% filter(Qtest_2$p_star>0.01)
ILS_0.01 <- clean_data(ILS_0.01)
ILS_0.01<- assign_high_median_low(ILS_0.01,
                                  cols = c("12|34","13|24","14|23"))


##################   针对每一个基因树进行所有三联体的提取   ######################


# 全局物种顺序
tnames <- taxonNames(sptree_list)
QT_sp <- get_qt_df(sptree_list,tnames)
QT_sp_2 <- QT_sp[,(ncol(QT_sp)-4):ncol(QT_sp)]

message("5. start loop")
# 树文件路径
tree_files <- list.files("../QuIBL_genetree_rt_final/", pattern = "\\.tre$", full.names = TRUE)
group_size <- 50
n_groups <- ceiling(length(tree_files) / group_size)

tree_groups <- list()
for(i in 1:n_groups) {
  start_idx <- (i - 1) * group_size + 1
  end_idx <- min(i * group_size, length(tree_files))
  tree_groups[[i]] <- tree_files[start_idx:end_idx]
}

tree_files_7 <- tree_groups[[7]]

#单个基因树
for(file in tree_files_7){
  tryCatch({
    print(file)
    
    gene_name <- tools::file_path_sans_ext(basename(file)) %>% sub("_.*", "",.)
    
    tr <- read.tree(file)
    tr <- list(tr)
    class(tr) <- "multiPhylo"
    
    # 当前基因树存在的物种
    present_taxa <- taxonNames(tr)
    
    # 若少于4个物种则跳过
    if (length(present_taxa) < 4) {
      warning(paste("Tree", basename(file), "has less than 4 taxa, skipped."))
      next
    }
    
    # 关键步骤：根据物种树顺序重新排列
    sorted_taxa <- tnames[tnames %in% present_taxa]
    
    # 调用 quartetTable
    df <- get_qt_df(tr, sorted_taxa)
    
    df_2 <- df %>% select(comb)
    
    #获得与物种树冲突的集合
    df_3 <- df_2 %>% 
      filter(!df_2$comb %in% QT_sp_2$comb)
    
    #从与物种树冲突的集合中再次挑选，不属于ILS也不属于IH的集合
    df_without <- df_3 %>%
      mutate(
        type = str_extract(comb, "12\\|34|13\\|24|14\\|23"),
        comb = str_trim(str_remove(comb, "12\\|34|13\\|24|14\\|23"))
      )%>% 
      filter(!comb %in% IH_0.01$comb & !comb %in% ILS_0.05$comb)
    
    df_without2 <- df_without %>%
      inner_join(all, by = "comb")
    
    write.csv(df_without2,paste0(gene_name,"_without_ILS_IH.csv"),row.names = F)
    
    #开始计算IH，ILS占比
    df_IH_0.01 <- get_contribution(df_3,IH_0.01)
    df_IH_0.001 <- get_contribution(df_3,IH_0.001)
    df_IH_0.0001 <- get_contribution(df_3,IH_0.0001)
    df_IH_0.00001 <- get_contribution(df_3,IH_0.00001)
    
    df_ILS_0.01 <- get_contribution(df_3,ILS_0.01)
    df_ILS_0.05 <- get_contribution(df_3,ILS_0.05)
    
    write.csv(df_IH_0.01,paste0(gene_name,"_df_IH_0.01.csv"),row.names = F)
    write.csv(df_IH_0.001,paste0(gene_name,"_df_IH_0.001.csv"),row.names = F)
    write.csv(df_IH_0.0001,paste0(gene_name,"_df_IH_0.0001.csv"),row.names = F)
    write.csv(df_IH_0.00001,paste0(gene_name,"_df_IH_0.00001.csv"),row.names = F)
    write.csv(df_ILS_0.01,paste0(gene_name,"_df_ILS_0.01.csv"),row.names = F)
    write.csv(df_ILS_0.05,paste0(gene_name,"_df_ILS_0.05.csv"),row.names = F)
    
    
    result <- data.frame(
      genename=gene_name,
      d_IH_0.01=nrow(df_IH_0.01)/nrow(df),
      d_IH_0.001=nrow(df_IH_0.001)/nrow(df),
      d_IH_0.0001=nrow(df_IH_0.0001)/nrow(df),
      d_IH_0.00001=nrow(df_IH_0.00001)/nrow(df),
      d_ILS_0.01=nrow(df_ILS_0.01)/nrow(df),
      d_ILS_0.05=nrow(df_ILS_0.05)/nrow(df),
      all_diff=nrow(df_3)/nrow(df),
      no_IH_0.01_ILS_0.05_but_diff=(nrow(df_3)-nrow(df_IH_0.01)-nrow(df_ILS_0.05))/nrow(df)
    )
    write.csv(result,paste0(gene_name,"_result.csv"),row.names = F)
  }, error = function(e) {
    warning(paste("Error processing", basename(file), ":", e$message))
  })
}