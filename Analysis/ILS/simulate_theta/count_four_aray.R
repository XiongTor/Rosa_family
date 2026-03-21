# # !/usr/bin/Rscript
# # Author: Tao Xiong
# # Date: 2025.05.04
# # Usage: This script is used to count the quartet frequencies from gene trees.

# # ==== Main Script Start ====

# #install.packages("MSCquartets")

library(ape)
library(Quartet)
library(dplyr)
library(ggplot2)
library(MSCquartets)
library(data.table)

files <- commandArgs(T)

#########################################  读树并修正  ##########################################
message("Reading species & gene trees...")

gtrees = read.tree(files[1])

sptree=read.tree(files[2])
sptree <- drop.tip(sptree,c("Elaeagnus_angustifolia","Zelkova_schneideriana","Morus_indica"))
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"

simulate_list = lapply(files[3:length(files)], read.tree)
names(simulate_list) = paste0("sim", seq_along(simulate_list))


# 提取亚科树，主要是作为参考，由于后续分离四联体需要multiPhylo格式树文件，因此在此处进行转换
Amy <- c("Prunus_padus","Exochorda_racemosa_subsp._serratifolia","Gillenia_trifoliata","Coleogyne_ramosissima","Kerria_japonica","Oemleria_cerasiformis","Prinsepia_uniflora","Rhodotypos_scandens","Lyonothamnus_floribundus","Alniaria_alnifolia","Amelanchier_laevis","Aria_edulis","Aronia_arbutifolia","Chaenomeles_speciosa","Chamaemeles_coriacea","Cormus_domestica","Cotoneaster_frigidus","Crataegus_hupehensis","Cydonia_oblonga","Dichotomanthes_tristaniicarpa","Eriobotrya_japonica","Griffitharia_hemsleyi","Hesperomeles_cuneata","Kageneckia_oblonga","Karpatiosorbus_bristoliensis","Lindleya_mespiloides","Macromeles_tschonoskii","Malacomeles_denticulata","Malus_domestica","Micromeles_alnifolia","Osteomeles_schweriniae","Peraphyllum_ramosissimum","Phippsiomeles","Photinia_prunifolia","Pourthiaea_amphidoxa","Pseudocydonia_sinensis","Pyracantha_coccinea","Pyrus_communis","Rhaphiolepis_ferruginea","Sorbus_aucuparia","Stranvaesia_nussia","Thomsonaria_caloneura","Torminalis_glaberrima","Vauquelinia_australis","Weniomeles_bodinieri","Neillia_sinensis","Physocarpus_opulifolius","Sorbaria_kirilowii_var._arborea","Aruncus_dioicus","Spiraea_blumei")
Dry <- c("Cercocarpus_ledifolius","Chamaebatia_foliolosa","Dryas_ajanensis_subsp._ajanensis","Purshia_tridentata")
Ros <- c("Acaena_ovalifolia","Agrimonia_pilosa_var._pilosa","Bencomia_exstipulata","Margyricarpus_pinnatus","Polylepis_tarapacana","Sanguisorba_minor","Coluria_longifolia","Geum_urbanum","Taihangia_rupestris","Waldsteinia_ternata","Alchemilla_faeroensis","Argentina_anserina","Chamaecallis_perpusilloides","Chamaerhodos_erecta","Comarum_palustre","Dasiphora_fruticosa","Drymocallis_arguta","Fragaria_nilgerrensis","Potaninia_mongolica","Potentilla_acaulis","Sibbaldia_parviflora","Sibbaldianthe_bifurca","Rosa_chinensis","Rubus_argutus","Filipendula_ulmaria")


subtree_Amy <- keep.tip(sptree, Amy)
subtree_Amy_list <- list(subtree_Amy)
class(subtree_Amy_list) <- "multiPhylo"

subtree_Dry <- keep.tip(sptree, Dry)
subtree_Dry_list <- list(subtree_Dry)
class(subtree_Dry_list) <- "multiPhylo"

subtree_Ros <- keep.tip(sptree, Ros)
subtree_Ros_list <- list(subtree_Ros)
class(subtree_Ros_list) <- "multiPhylo"


# 修正模拟树末端枝和负枝长
simulate_list <- lapply(simulate_list, function(tree_set) {
  lapply(tree_set, function(tr) {
    Ntip <- length(tr$tip.label)
    # ≤0的枝长设为极小值
    tr$edge.length[tr$edge.length <= 0] <- 1e-4
    return(tr)
  })
})

############################   获得四连体组合数据框  ###############################
### 得到数据集中所有物种名字
tnames=taxonNames(sptree_list)

### 函数：生成 QT_df + comb列
get_qt_df <- function(tree,names){
  QT = quartetTable(tree,names, progressbar = TRUE)
  df = as.data.frame(QT)
  df$comb <- apply(df, 1, function(x){
    taxa <- names(x)[which(x == 1)]
    paste(taxa, collapse = " ")
  })
  return(df)
}

#提前准备需要读出的名字，用于保存进度，如不需要可以注释掉
qt_cache_file <- "QT_cache.RDS"

message("Calculating quartet tables...")

#准备参考树及参考树中四联体组合
QT_sp <- get_qt_df(sptree_list,tnames)

tname_Amy <- taxonNames(subtree_Amy_list)
QT_Amy <- get_qt_df(subtree_Amy_list,tname_Amy)

tname_Dry <- taxonNames(subtree_Dry_list)
QT_Dry <- get_qt_df(subtree_Dry_list,tname_Dry)

tname_Ros <- taxonNames(subtree_Ros_list)
QT_Ros <- get_qt_df(subtree_Ros_list,tname_Ros)

# #基因树集合及四连体运算
QT_real   <- get_qt_df(gtrees,tnames)

#模拟树集合及四连体运算
QT_sim_df_list <- lapply(simulate_list, function(tree_set) {
  get_qt_df(tree_set, tnames)
})

# 给 QT_sim_df_list 列加后缀，如“_1”，"_2"
for(nm in names(QT_sim_df_list)){
   colnames(QT_sim_df_list[[nm]])[1:3] <- paste0(colnames(QT_sim_df_list[[nm]])[1:3], "_", nm)
}

# 保存缓存
message("QT data saved to QT_cache.RDS")
saveRDS(list(QT_real=QT_real, QT_sim_df_list=QT_sim_df_list), qt_cache_file)

#读取之前保存的文件
# qt_cache <- readRDS("QT_cache.RDS")
# QT_real <- qt_cache$QT_real

   
#######################   提取所需数据并计算占比   #############################
#处理作为参考的物种树，提取需要的列，并更改comb列的内容，使得程度可用
QT_sp_2 <- QT_sp[,(ncol(QT_sp)-4):ncol(QT_sp)]
QT_sp_2$comb <- sapply(strsplit(QT_sp_2$comb, "\\s+"), function(x) {
  paste(x[1:4], collapse = " ")
})

QT_Amy_2 <- QT_Amy[,(ncol(QT_Amy)-4):ncol(QT_Amy)]
QT_Amy_2$comb <- sapply(strsplit(QT_Amy_2$comb, "\\s+"), function(x) {
  paste(x[1:4], collapse = " ")
})

QT_Dry_2 <- QT_Dry[,(ncol(QT_Dry)-4):ncol(QT_Dry)]
QT_Dry_2$comb <- sapply(strsplit(QT_Dry_2$comb, "\\s+"), function(x) {
  paste(x[1:4], collapse = " ")
})

QT_Ros_2 <- QT_Ros[,(ncol(QT_Ros)-4):ncol(QT_Ros)]
QT_Ros_2$comb <- sapply(strsplit(QT_Ros_2$comb, "\\s+"), function(x) {
  paste(x[1:4], collapse = " ")
})

#将所有参考树的结果保存到一个向量中，方便后续选择
reference <- c("QT_sp_2", "QT_Amy_2", "QT_Ros_2")

#处理基因树集合，提取需要的列，并更改comb列的内容，使得程度可用
QT_real_2 <- QT_real[,(ncol(QT_real)-4):ncol(QT_real)]
QT_real_2$comb <- sapply(strsplit(QT_real_2$comb, "\\s+"), function(x) {
  paste(x[1:4], collapse = " ")
})


#定义函数，以参考列为基准，提取出其它数据框中对应四连体的比例
get_prop_compare_fast <- function(ref_df, target_df){
  qcols <- c("12|34","13|24","14|23","1234")
  dt <- data.table(merge(ref_df, target_df, by="comb", suffixes=c("_ref","_tgt")))
  # 1. 参考数据框行中哪列为1
  ref_mat <- as.matrix(dt[, paste0(qcols, "_ref"), with=FALSE])
  dt[, type := qcols[max.col(ref_mat, ties.method="first")]]
  # 2. 目标数据总和
  tgt_mat <- as.matrix(dt[, paste0(qcols, "_tgt"), with=FALSE])
  total <- rowSums(tgt_mat)
  # 3. proportion = 目标data中对应topology / 总和
  dt[, prop := tgt_mat[cbind(1:.N, match(type, qcols))] / total]
  # 保留comb、原12|34等列、type、prop
  keep_cols <- c(
    "comb",
    paste0(qcols, "_ref"),
    paste0(qcols, "_tgt"),
    "type",
    "prop"
  )
  return(dt[, ..keep_cols])
}

for(ref_name in reference){
  ref <- get(ref_name)
  cat("Processing with reference:", ref_name, "\n")
#开始以参考树为基准，处理基因树集合
result_genetree <- get_prop_compare_fast(ref, QT_real_2)

#开始以参考树为基准，处理模拟树集合
# QT_sim_df_list <- qt_cache$QT_sim_df_list

QT_sim_df_list_sub <- lapply(names(QT_sim_df_list), function(nm){
  # 提取子集，并更改comb列的内容
  df <- QT_sim_df_list[[nm]]
  mm <- df[, (ncol(df)-4):ncol(df)]
  mm$comb <- sapply(strsplit(mm$comb, "\\s+"), function(x) {
    paste(x[1:4], collapse = " ")
  })
  #开始以参考树为基准，处理基因树集合
  result <- get_prop_compare_fast(ref, mm)
  # 只给非 comb 的列加后缀
  colnames(result) <- ifelse(colnames(result) == "comb",
                         "comb",
                         paste0(colnames(result), "_", nm))
  result
})

#给list中对应的sim重命名加上后缀编号
names(QT_sim_df_list_sub) <- names(QT_sim_df_list)


#################################  计算相关性   #######################################

# 对真实数据，从基因树集合提取所需列
real_df <- result_genetree[, c("comb", "prop")]

# 创建保存结果的列表
lm_results <- list()
plot_list <- list()

# 遍历模拟列表,并与真实数据进行线性回归集合判断相关性，同时作图展示结果
 for(nm in names(QT_sim_df_list_sub)){
  
  cat("处理:", nm, "\n")
  
  # 提取模拟 df数据框
  sim_df <- as.data.frame(QT_sim_df_list_sub[[nm]])
  
  # 找到 prop 列名，例如 prop_sim1 / prop_sim2 ...
  prop_col <- grep("^prop_", colnames(sim_df), value = TRUE)
  
  if(length(prop_col) != 1){
    stop(paste("在", nm, "中找不到或找到多于1个prop列"))
  }
 
  # 按comb列进行merge
  df <- merge(real_df, sim_df[, c("comb", prop_col)], by = "comb")
  colnames(df)[colnames(df) == prop_col] <- "prop_sim"
  
  # 线性回归
  fit <- lm(prop_sim ~ prop, data = df)
  lm_results[[nm]] <- summary(fit)
  
  # 散点抽样
  n_points <- min(2000, nrow(df))   # 最多 2000 个点
  df_sample <- df[sample(nrow(df), n_points), ]
  
  # 画图
  p <- ggplot(df_sample, aes(x = prop, y = prop_sim)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], color="red") + 
    theme_bw() +
    labs(title = paste("Real vs", nm),
         y = paste0("Simulated: ", nm),
         x = "Real prop",
         subtitle = paste("R² =", round(summary(fit)$r.squared, 4)))
  
  plot_list[[nm]] <- p
  
  # 保存图
  ggsave(filename = paste0(ref_name,"_scatter_", nm, ".pdf"), plot = p, width = 6, height = 9)
 }
}

# 输出回归结果到具体表格中
# 定义输出文件
# outfile <- "lm3_results.txt"
# # 打开连接
# sink(outfile) 
# for(nm in names(lm_results)){
#   cat("\n===== ", nm, " =====\n")
#   print(lm_results[[nm]])
# }
# sink()