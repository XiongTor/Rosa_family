#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-07-22
# Description: Plot densitree to display the gene trees distinction

# ==== Main Script Start ====

library(devtools)
# source("https://bioconductor.org/biocLite.R")
library(tidyr)
library(dplyr)
library(phytools)
library(phangorn)
library(ggtree)
library(ape)
library(bindrcpp)
library(scales)

# Load the tree
sptree <- read.tree("data_tree/rosa_tribe_Astral_species_rt_e_nolength.tre")
genetree <- read.tree("data_tree/rosa_genetrees_e_nolength.tre")

is.binary(sptree)


## 过滤
# 统计所有基因树的枝长，去掉其中枝长太长的树
branch_lengths <- sapply(genetree, function(tr) sum(tr$edge.length))

# 用 IQR 方法定义“异常高”的阈值
Q1 <- quantile(branch_lengths, 0.25)
Q3 <- quantile(branch_lengths, 0.75)
IQR_val <- Q3 - Q1
threshold <- Q3 + 1 * IQR_val

# 找出不是异常值的索引
keep_idx <- which(branch_lengths <= threshold)

# 筛选出保留的树
filtered_trees <- genetree[keep_idx]
class(filtered_trees) <- "multiPhylo"

# 输出删除了多少棵树
cat("删除了", length(genetree) - length(filtered_trees), "棵异常长的树\n")

# write.tree(filtered_trees,"data_tree/rosa_genetrees_e_clean.tre")
# 清理所有树的tips，只保留物种树中的tips 
cleaned_trees <- lapply(filtered_trees, function(tr) {
  common_tips <- intersect(tr$tip.label, sptree$tip.label)
  drop.tip(tr, setdiff(tr$tip.label, common_tips))
})
class(cleaned_trees) <- "multiPhylo"


# 如果树太多可以抽取其中一部分用于建树
if (length(cleaned_trees) > 353) {
  cleaned_trees_2 <- cleaned_trees[sample(1:length(cleaned_trees), 353)]
}
class(cleaned_trees) <- "multiPhylo"




## 排列物种顺序
desired_order <- rev(c("Pyrus_communis","Malus_domestica","Sorbus_aucuparia","Gillenia_trifoliata","Lyonothamnus_floribundus","Prunus_padus","Rhodotypos_scandens","Oemleria_cerasiformis","Spiraea_purpurea","Physocarpus_opulifolius","Dryas_drummondii","Potentilla_acaulis","Agrimonia_pilosa","Rosa_chinensis","Geum_urbanum","Rubus_idaeus","Filipendula_ulmaria","Morus_alba","Zelkova_schneideriana","Elaeagnus_pungens_hangzhou"))

# 对共识树进行重新排列
sptree <- rotateConstr(sptree, desired_order)


# ### group
# subf <- read.csv("Rosaceae_genus_accepted_lastest.csv")
# tips <- as.data.frame(sptree$tip.label)
# colnames(tips) <- "species"
# tips$genus <- sapply(strsplit(tips$species,"_"),`[`,1)
# 
# groups <- merge(tips,subf,by.x = "genus",by.y = "Genera",all.x = T) %>% 
#   select(genus,species,Subfamilies)
# groups$Subfamilies[is.na(groups$Subfamilies)] <- "outgroups"
# 
# #分配颜色
# tip_colors <- sapply(tips$species, function(tip) {
#   if (tip %in% groups$species[groups$Subfamilies == "Rosoideae"]) {
#     "blue"
#   } else if (tip %in% groups$species[groups$Subfamilies == "Amygdaloideae"]) {
#     "red"
#   } else if (tip %in% groups$species[groups$Subfamilies == "Dryadoideae"]) {
#     "yellow"
#   } else {
#     "black"
#   }
# })
# names(tip_colors) <- tips$species
# 
# colors_list <- lapply(cleaned_trees, function(tree) {
#   unname(tip_colors[tree$tip.label])
# })
# 
# colors_list <- as.vector(colors_list)
### 重复物种树多次，然后读入
spmore <- read.tree("data_tree/rosa_tribe_Astral_species_rt_e_nolength_cpmore.tre")

# 绘制 DensiTree 图

densiTree(cleaned_trees,
          consensus = sptree,
          type = "cladogram",
          label.offset=2,
          width=1,
          alpha = 0.08,
          scaleX =T,
          direction = "rightwards",
          col = "#272361",
          scale.bar = F)
par(new=TRUE)
densiTree(spmore,
          consensus = sptree,
          type = "cladogram",
          label.offset=0.01,
          width=3.5,
          alpha = 1,
          scaleX =T,
          direction = "rightwards",
          col = "black",
          scale.bar = F)




       