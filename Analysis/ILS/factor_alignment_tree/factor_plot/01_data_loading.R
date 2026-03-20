# !/usr/bin/Rscript
# Author: XiongTao
# Description: Data loading and preprocessing

library(dplyr)
library(tidyverse)

# Load chloroplast data
file_c <- list.files("./chloroplast/")

df_chloroplast <- list()
for (name in file_c) {
  df_chloroplast[[name]] <- read.table(paste0("./chloroplast/", name), header = T)
  df_chloroplast[[name]]$Alignment_name <- sapply(
    strsplit(df_chloroplast[[name]]$Alignment_name, "\\."), `[`, 1
  )
}

merged_chloroplast <- reduce(df_chloroplast, full_join, by = "Alignment_name") %>%
  select(Alignment_name, No_of_taxa, Alignment_length, GC_content,
         Proportion_variable_sites, Proportion_parsimony_informative,
         Mean_evolutionary_rate, Mean_evolutionary_rate_without0, RCV,
         Mean_bipartition, Mean_internal_branch, Mean_teminal_branch,
         treeness, saturation, absolute_saturation, treeness.RCV)

# Load nuclear gene data
# 如果需要对数据进行挑选，则将挑选后的基因树放置在对应文件夹下后进行读取，利用新读取到的名称进行筛选
tt <- list.files(path = "./tree/orthofinder_genes/",pattern = "\\.tre",full.names = T)

# 提取路径
Alignment_name <- basename(tt)
df <- data.frame(Alignment_name)
df$Alignment_name <- sub("\\_rt.tre$","",df$Alignment_name)

# 读取所有的数据并进行后续处理
file_path <- list.files(path = "./orthofinder/",pattern = "\\.txt",full.names = TRUE)

all_data <- lapply(file_path, read.table, header = TRUE, sep = "\t")

data_list <- c(list(df),all_data)

final <- Reduce(function(x,y) merge(x,y,by="Alignment_name",all.x=T),data_list)

merged_orthofinder <- final %>%
  select(Alignment_name, No_of_taxa, Alignment_length, GC_content,
         Proportion_variable_sites, Proportion_parsimony_informative,
         Mean_evolutionary_rate, Mean_evolutionary_rate_without0, RCV,
         Mean_bipartition, Mean_internal_branch, Mean_teminal_branch,
         treeness, saturation, absolute_saturation, treeness.RCV)

# Add dataset labels
merged_orthofinder$dataset <- "orthofinder"
merged_chloroplast$dataset <- "chloroplast"

# Save processed data
saveRDS(merged_chloroplast, "merged_chloroplast.rds")
saveRDS(merged_orthofinder, "merged_orthofinder.rds")

cat("Data loading completed.\n")
cat("Chloroplast genes:", nrow(merged_chloroplast), "\n")
cat("Nuclear genes:", nrow(merged_orthofinder), "\n")
