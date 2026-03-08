#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2026-02-08
# Description: plot for quibl

# ==== Main Script Start ====


# install.packages("devtools")
# library(devtools)
# install_github("nbedelman/quiblR")
# library(quiblR)

library("quiblR")
library("ggplot2")
library("ape")
library("hash")
library("ggtree")
library("ggpubr")
library("dplyr")
library("patchwork")

#加载数据
speciesTree  <- read.tree("rosa_tribe_Astral_species.rt.tre")
quiblOutput  <- read_csv_quibl("results_qulbl.txt")
originalTrees  <- read.tree("all_4_trees.tre")
# originalTrees <- phytools::read.newick("rosa_ags353_4_trees.tre")

#判断拓扑是否与物种树一致，1表示不一致，0表示一致
message("1. Judge whether the quartet tree topology is the same as the species tree.")
quiblOutput <- mutate(quiblOutput, isDiscordant=as.integer(!apply(quiblOutput, 1, isSpeciesTree, sTree=speciesTree)))

#计算每个triplet的渗入概率的总和
totalTrees <- sum(quiblOutput$count)/length(unique(quiblOutput$triplet))
quiblOutput <- mutate(quiblOutput,totalIntrogressionFraction=(mixprop2*count*isDiscordant)/totalTrees)

#汇总，同时新加入一列判断是否具有显著差异，应该是以ΔBIC值是否>10来作为判断条件，判断渐渗+ILS模型是否显著优于only_ILS模型。、
#BIC越小数据集合越好
largeQuiblOutput <- mutate(quiblOutput,
                           isDiscordant = as.integer(! apply(quiblOutput, 1, isSpeciesTree, sTree=speciesTree)),
                           isSignificant = as.integer(apply(quiblOutput, 1, testSignificance, threshold=10)),
                           totalIntrogressionFraction=(mixprop2*count*isDiscordant)/totalTrees)


#################################### getIntrogressionSummary 函数详解  ####################################
# 如果所有基因树拥有的物种数量是一致的，那么可以不用重新定义这个函数，直接使用原R包定义的此函数即可，其余情况建议还是重新编译一下此函数再使用。

#用于测试函数
quiblOut <- quiblOutput

getIntrogressionSummary2 <- function (quiblOut, speciesTree, summaryType = "mean") 
{
  species <- unique(quiblOut[, "outgroup"])
  orderedSpecies <- factor(species, levels = speciesTree$tip.label)
  correlationDict <- makeEmptyHash(orderedSpecies)
for (row in 1:nrow(quiblOut)){
  #计算第一种三联体中的基因树总数量
  numTrees <- sum(quiblOut[which(quiblOut$triplet==quiblOut$triplet[row]),]$count)
  #提取出除外类群外的其它物种
  taxa <- setdiff(unlist(strsplit(as.character(quiblOut[row,][,"triplet"]),"_")),quiblOut[row,][,"outgroup"])
  #制造两个物种的排列组合
  key1 <- paste(sort(taxa)[1], sort(taxa)[2],sep = "_")
  key2 <- paste(sort(taxa)[2], sort(taxa)[1],sep = "_")
  #判断是否与物种树一致，从而来判断是否需要计算渐渗概率
  if(! isSpeciesTree(quiblOut[row,],speciesTree)){
    #判断数据到底加到哈希表中的哪个组合中
    if (key1 %in% keys(correlationDict)){
      correlationDict[[key1]] <- append(correlationDict[[key1]],quiblOut[row,][,"mixprop2"]*quiblOut[row,][,"count"]/numTrees)
    } else {
      correlationDict[[key2]] <- append(correlationDict[[key2]],quiblOut[row,][,"mixprop2"]*quiblOut[row,][,"count"]/numTrees)
    }
  }
}

inp <- data.frame(tax1=character(),tax2=character(), value=numeric())
#计算两个物种对之间的平均值，并重复两次方便画热图
for(k in keys(correlationDict)){
  taxa <- unlist(strsplit(k,"_"))[1:2]
  inp <- rbind(inp, data.frame(tax1=taxa,
                               tax2=rev(taxa),
                               value=rep(mean(correlationDict[[k]]),2)))
}
#整理最终表格的结果与系统发育树一致
inp[,"tax1"] <- factor(inp[,"tax1"], levels = speciesTree$tip.label)
inp[,"tax2"] <- factor(inp[,"tax2"], levels = speciesTree$tip.label)
#把na值填充为0
inp[,"value"][is.na(inp[,"value"])] <- 0
#前面已经转化为因子了，这里是按照系统发育树进行排序
inp <- inp[order(inp[["tax1"]], inp[["tax2"]]),]
return(inp)
}


########################################  按照如上逻辑，计算ILS的概率  ##################################
getILSSummary <- function(quiblOut,speciesTree,summaryType="mean"){
  species <- unique(quiblOut[,"outgroup"])
  orderedSpecies <- factor(species, levels=speciesTree$tip.label)
  correlationDict <- makeEmptyHash(orderedSpecies)
  
  for (row in 1:nrow(quiblOut)){
    numTrees <- sum(quiblOut[which(quiblOut$triplet==quiblOut$triplet[row]),]$count)
    taxa <- setdiff(unlist(strsplit(as.character(quiblOut[row,][,"triplet"]),"_")),quiblOut[row,][,"outgroup"])
    key1 <- paste(sort(taxa)[1], sort(taxa)[2],sep = "_")
    key2 <- paste(sort(taxa)[2], sort(taxa)[1],sep = "_")
    if(! isSpeciesTree(quiblOut[row,],speciesTree)){
      if (key1 %in% keys(correlationDict)){
        correlationDict[[key1]] <- append(correlationDict[[key1]],quiblOut[row,][,"mixprop1"]*quiblOut[row,][,"count"]/numTrees)
      } else {
        correlationDict[[key2]] <- append(correlationDict[[key2]],quiblOut[row,][,"mixprop1"]*quiblOut[row,][,"count"]/numTrees)
      }
    }
  }
  
  inp <- data.frame(tax1=character(),tax2=character(), value=numeric())
  for(k in keys(correlationDict)){
    taxa <- unlist(strsplit(k,"_"))[1:2]
    inp <- rbind(inp, data.frame(tax1=taxa,
                                 tax2=rev(taxa),
                                 value=rep(mean(correlationDict[[k]]),2)))
  }
  inp[,"tax1"] <- factor(inp[,"tax1"], levels = speciesTree$tip.label)
  inp[,"tax2"] <- factor(inp[,"tax2"], levels = speciesTree$tip.label)
  inp[,"value"][is.na(inp[,"value"])] <- 0
  inp <- inp[order(inp[["tax1"]], inp[["tax2"]]),]
  return(inp)
}

############################### 渐渗概率画热图 #################################

#########################  画出树图并提取其中的tips排列顺序 #################

#### 画出树图，删去哪些在热图中没有的tips，并使得tips排列方式与热图一致
message("2. Plot species tree")
speciesTree$edge.length <- speciesTree$edge.length / 10
speciesTree <- ladderize(speciesTree)
taxa_order <- unique(introgressionSummary$tax2)

speciesTreeSubset <- ggtree(extractTripletTree(speciesTree, unique(introgressionSummary$tax2)))+
  geom_tiplab(align = TRUE,size=2.5)+
  coord_cartesian(clip = "off") + 
  theme(plot.margin = margin(10, 100, 10, 10))

#按照绘图的tips顺序提取树中的tips名字
tip_order <- speciesTreeSubset$data %>%
  filter(isTip == TRUE) %>%
  arrange(desc(y)) %>% # 按y坐标降序排列
  pull(label)

#明确最终热图中使用的tips的顺序
taxa_in_heatmap <- unique(c(as.character(introgressionSummary$tax1), as.character(introgressionSummary$tax2)))
final_taxa_order <- tip_order[tip_order %in% taxa_in_heatmap]

####################### 结束 #################################################


################################  渐渗数据集合 #############################
#计算与物种树不一致的基因树中，物种对之间渐渗概率的拟合，加权平均数
message("3. Calculate the probability of integration.")
introgressionSummary <- getIntrogressionSummary2(quiblOutput,speciesTree)

message("3.1 Save data -- integration.")
# saveRDS(introgressionSummary2,"introgressionSummary2.RDS")
# tt <- readRDS("introgressionSummary2.RDS")

# 创建一个新的数据框副本
introgressionSummary_ordered <- introgressionSummary

# 将 tax1 和 tax2 列转换为因子，并使用从树图得到的顺序 (final_taxa_order)
# 这会强制 ggplot 按照你指定的顺序来绘制坐标轴
introgressionSummary_ordered$tax1 <- factor(introgressionSummary_ordered$tax1, levels = rev(final_taxa_order))

# 对于y轴，顺序是反的，所以我们需要使用 rev() 来反转顺序
introgressionSummary_ordered$tax2 <- factor(introgressionSummary_ordered$tax2, levels = rev(final_taxa_order))



################################ 以同样的方法调整ILS数据集合 ######################
message("4. Calculate the probability of ILS.")
ILS <- getILSSummary(quiblOutput,speciesTree)

message("4.1 Save data -- integration.")
saveRDS(ILS,"ILS.RDS")
#tt <- readRDS("ILS.RDS")

# 创建一个新的数据框副本
ILS_ordered <- ILS

# 将 tax1 和 tax2 列转换为因子，并使用从树图得到的顺序 (final_taxa_order)
# 这会强制 ggplot 按照你指定的顺序来绘制坐标轴
ILS_ordered$tax1 <- factor(ILS_ordered$tax1, levels = rev(final_taxa_order))

# 对于y轴，顺序是反的，所以我们需要使用 rev() 来反转顺序
ILS_ordered$tax2 <- factor(ILS_ordered$tax2, levels = rev(final_taxa_order))



################################ 画图 ##############################
message("5. plot hotmap")

# 设置索引列用于讲最终热图左上画为ILS右下画为IH
tax1_index_ils <- match(as.character(ILS_ordered$tax1), rev(final_taxa_order))
tax2_index_ils <- match(as.character(ILS_ordered$tax2), rev(final_taxa_order))

tax1_index_intro <- match(as.character(introgressionSummary_ordered$tax1), rev(final_taxa_order))
tax2_index_intro <- match(as.character(introgressionSummary_ordered$tax2), rev(final_taxa_order))

# 添加辅助列
ILS_ordered$region <- ifelse(tax1_index_ils < tax2_index_ils, "upper", "ignore")
introgressionSummary_ordered$region <- ifelse(tax1_index_intro > tax2_index_intro, "lower", "ignore")


#根据计算结果画出热图
summaryMatrix <- ggplot() +
  # 左上角画 ILS
  geom_tile(data = ILS_ordered,
            aes(x = tax1, y = tax2,
                fill = ifelse(region == "upper", value,  NA_real_)
                ),
            color = "white") +
  #右下角画 introgression
  geom_tile(data = introgressionSummary_ordered[introgressionSummary_ordered$region == "lower", ],
            aes(x = tax1, y = tax2, fill = value),
            color = "white") +
  scale_fill_gradient2(low = "white", high = "red", mid = "#FFD700", na.value = "grey100",
                       midpoint = max(ILS_ordered$value)/2, limit = c(0, max(ILS_ordered$value)),
                       name="Average ILS/IH fraction") +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = seq(1.5, nrow(ILS_ordered) + 0.5, 1), alpha = 0.7) +
  geom_hline(yintercept = seq(1.5, nrow(ILS_ordered) + 0.5, 1), alpha = 0.7) +
  labs(x = "", y = "") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90,hjust=1,vjust = 0.5,size=7.5)
  )+
  theme(plot.margin = margin(10, 10, 10, 0))


#汇总

ggarrange(speciesTreeSubset, summaryMatrix,
            ncol = 2,nrow=1,widths=c(1,2.4),align = "h",common.legend = T)

pdf("summary_heatMap.pdf",width = 10, height = 8,family = "Helvetica")       # 避免中文乱码；如需中文字体改成 "GB1"
ggarrange(speciesTreeSubset, summaryMatrix,
          ncol = 2, nrow = 1,
          widths = c(1, 2.4),
          align = "h",
          common.legend = TRUE)
dev.off()


############################ 画ILS/渐渗分布的密度图 ############################
## 自定义函数，对原本的基因树集合进行筛选，选择具有特定三联体的基因树

filterTreesByTriplet <- function(treeList, triplet) {
  # 拆分 triplet 名字
  triplet_species <- unlist(strsplit(triplet, "_"))
  
  # 保留包含所有3个物种的树
  filtered <- treeList[sapply(treeList, function(tree) {
    all(triplet_species %in% tree$tip.label)
  })]
  
  return(filtered)
}

#给quiblOutput添加两列，用于计算BIC值的大小
quiblOutput <- mutate(quiblOutput,
  ILS_BIC = BIC1Dist-BIC2Dist,
  IH_BIC = BIC2Dist-BIC1Dist
)

# 待检测的一对物种和所有基因树的外类群，这两种都是需要被去除的类型

parts <- c("Pyrus.communis", "Sorbus.aucuparia","Zelkova.schneideriana")

# 其它物种
others <- setdiff(tip_order,parts)





#穷尽所有组合
plot_list <- list()  # 存储所有图
index <- 1   
for(name in others){
  targets <- c("Pyrus.communis", "Sorbus.aucuparia", name)

#匹配到对应的三联体
matches <- sapply(strsplit(as.character(quiblOutput$triplet), "_"), function(taxa) {
  all(targets %in% taxa)
})

#挑出BIC插值更小的行========================
target_triplet <- quiblOutput[matches, ] %>% filter(isDiscordant=="1",!is.na(mixprop2), mixprop2 != 0) %>% slice(which.min(ILS_BIC))

trip <- as.character(target_triplet$triplet)

#过滤出含有这个triplet的基因树
filteredTrees <- filterTreesByTriplet(originalTrees, trip)

#得到branchLength的概率分布
LocusStats <- getPerLocusStats(quiblOutput = quiblOutput, trip = trip, treeList = filteredTrees, overallOut="Zelkova.schneideriana")

#计算去除部分极值后的最大值是多少
threshold <- quantile(LocusStats$branchLength, 0.99, na.rm = TRUE)
max_trimmed <- max(LocusStats$branchLength[LocusStats$branchLength <= threshold], na.rm = TRUE)

# max_trimmed <- max(LocusStats$branchLength)



#集合不同模型的曲线
ILSOnly <- getILSOnlyDist(0,max_trimmed,subset(quiblOutput, outgroup == as.character(target_triplet$outgroup) & triplet==trip))
ILSMix <- getILSMixtureDist(0,max_trimmed,subset(quiblOutput, outgroup == as.character(target_triplet$outgroup) & triplet==trip))
nonILSMix <- getNonILSMixtureDist(0,max_trimmed,subset(quiblOutput, outgroup==as.character(target_triplet$outgroup) & triplet==trip))

#画图
P_ILS <- ggplot() +
  #第一层：直方图
  geom_histogram(data = LocusStats,
                 aes(x = branchLength),
                 bins = round(30/0.2*max_trimmed,0),
                 fill = "gray50") +
  # 第二层：拟合曲线
  geom_line(data = ILSOnly,
            aes(x = x, y = y),
            color = "red", size = 1.5, alpha=0.7,linetype = "dashed") +
  labs(x = "branchLength", y = "count",title=name)+
  xlim(0, max_trimmed) +
  theme_bw()+
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )+
  annotate("text", 
           x = max_trimmed,  # x轴位置在最大处
           y = Inf,  # y轴设为Inf表示顶部=========================
           label = paste0("Only ILS            \nILS_df_BIC = ", round(target_triplet$ILS_BIC, 2)), 
           hjust = 1.1,  # 水平方向稍微往左偏
           vjust = 1.5,  # 垂直方向往下偏一点
           size = 6)

P_IH_ILS <- ggplot() +
  #第一层：直方图
  geom_histogram(data = LocusStats,
                 aes(x = branchLength),
                 bins = round(30/0.2*max_trimmed,0),
                 fill = "gray50") +
  # 第二层：拟合曲线
  geom_line(data = ILSMix,
            aes(x = x, y = y),
            color = "red", size = 1.5, alpha=0.7,linetype = "dashed") +
  geom_line(data = nonILSMix,
            aes(x = x, y = y),
            color = "blue", size = 1.5, alpha=0.7,linetype = "dashed") +
  labs(x = "branchLength", y = "count",title=name)+
  xlim(0, max_trimmed) +
  theme_bw()+
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )+
  annotate("text", 
           x = max_trimmed,  # x轴位置在最大处
           y = Inf,  # y轴设为Inf表示顶部===========================
           label = paste0("IH+ILS          \nILS_df_BIC = ", round(target_triplet$ILS_BIC, 2)), 
           hjust = 1.1,  # 水平方向稍微往左偏
           vjust = 1.5,  # 垂直方向往下偏一点
           size = 6)

plot_list[[index]] <- P_ILS
index <- index + 1
plot_list[[index]] <- P_IH_ILS
index <- index + 1
}


plots_per_page <- 12  # 每页 3 行 × 4 列
n_pages <- ceiling(length(plot_list) / plots_per_page)

pdf("ILS_Pyrus.communis_Sorbus.aucuparia_251102.pdf", width = 16, height = 12)  # 打开PDF设备

for (i in seq_len(n_pages)) {
  start <- (i - 1) * plots_per_page + 1
  end <- min(i * plots_per_page, length(plot_list))
  
  p_page <- wrap_plots(plot_list[start:end], nrow = 3, ncol = 4, byrow = TRUE)+
    plot_annotation(title = "Pyrus.communis_Sorbus.aucuparia")
  print(p_page)  # 输出当前页图像
}

dev.off()


pdf("test.pdf", width = 16, height = 12)
marrangeGrob(grobs = plot_list, nrow = 3, ncol = 4) |> print()
dev.off()








LocusStats <- getPerLocusStats(quiblOutput = quiblOutput, trip = "Malus.domestica_Pyrus.communis_Sorbus.aucuparia", treeList = filteredTrees, overallOut="Zelkova.schneideriana")
head(LocusStats)

Hhsa_Htel_Hera_ILSOnly <- getILSOnlyDist(0,0.5,subset(quiblOutput, outgroup=="Malus.domestica" & triplet=="Malus.domestica_Pyrus.communis_Sorbus.aucuparia"))
Hhsa_Htel_Hera_ILSMix <- getILSMixtureDist(0,0.5,subset(quiblOutput, outgroup=="Malus.domestica" & triplet=="Malus.domestica_Pyrus.communis_Sorbus.aucuparia"))
Hhsa_Htel_Hera_nonILSMix <- getNonILSMixtureDist(0,0.5,subset(quiblOutput, outgroup=="Malus.domestica" & triplet=="Malus.domestica_Pyrus.communis_Sorbus.aucuparia"))


ggplot() +
  #第一层：直方图
  geom_histogram(data = LocusStats,
                 aes(x = branchLength),
                 bins = 210,
                 fill = "gray50") +
  # 第二层：拟合曲线
  geom_line(data = Hhsa_Htel_Hera_ILSOnly,
            aes(x = x, y = y),
            color = "red", size = 1.5, alpha=0.7,linetype = "dashed") +
  # geom_line(data = Hhsa_Htel_Hera_nonILSMix,
  #           aes(x = x, y = y),
  #           color = "blue", size = 1.5, alpha=0.7,linetype = "dashed") +
  labs(x = "branchLength", y = "count")+
  xlim(0, 0.5) +
  theme_bw()+
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )+
  annotate("text", 
           x = 0.5,  # x轴位置在最大处
           y = Inf,  # y轴设为Inf表示顶部
           label = "ΔBIC = -11", 
           hjust = 1.1,  # 水平方向稍微往左偏
           vjust = 1.5,  # 垂直方向往下偏一点
           size = 10)