# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2026.01.28
# This script is used to evaluate the contributions of ILS and IH.

library(ape)
library(phangorn)
library(MSCquartets)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(doParallel)
library(tidyr)
library(reshape2)
library(ggplot2)


# 0. 设置基因树的路径
# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1) {
  stop("Usage: Rscript script.R <tree_directory> [output_directory]
  Example: Rscript script.R ./genetrees_batch1/ ./output/")
}

tree_dir <- args[1]
output_dir <- ifelse(length(args) >= 2, args[2], "./12-MSCquartet/ags353/get_gene_contribution/")

# 创建输出目录
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

message(paste("Input directory:", tree_dir))
message(paste("Output directory:", output_dir))


# 1. 读入物种树 ----------------------------------------------------------------
message("1. read tree.")
sptree=read.tree("../rosa_ags353_treeshrink_sp_rt_oneoutg_final.tre")
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"
tnames <- taxonNames(sptree_list)


# 2. 读入MSCquartet的运算结果 --------------------------------------------------
message("2. load data")

load("../quartet_analysis_results.RData")

## 筛选出在不同限制条件下的渐渗物种和高ILS物种
Qtest_2 <- as.data.frame(Qtest_2)

# 3. 自定义函数 ----------------------------------------------------------------
message("3. load function")

## 为单个树计算QT矩阵，并生成comb列
get_qt_df <- function(tree, names,exclude_cols = c("12|34","13|24","14|23")){
  QT = quartetTable(tree, names, progressbar = TRUE)
  df = as.data.frame(QT)
  df$comb <- apply(df, 1, function(x){
    all_colnames <- names(x)
    target_col <- which(!all_colnames %in% exclude_cols)
    taxa <- all_colnames[target_col][which(x[target_col]==1)]
    paste(taxa, collapse = " ")
  })
  df$comb_all <- apply(df, 1, function(x){
    taxa <- names(x)[which(x == 1)]
    paste(taxa, collapse = " ")
  })
  df_2 <- df[,(ncol(df)-5):ncol(df)]
  return(df_2)
}


## 对已有的QT矩阵进行清理，新增comb列，同时去除其它无关信息
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



## 函数，用于判断大小
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


# 定义最终画图矩阵中的物种排列顺序
# 定义物种顺序
species_order <- c(
  "Outgroup",
  "Sibbaldia_parviflora",
  "Sibbaldianthe_bifurca",
  "Alchemilla_faeroensis",
  "Comarum_palustre",
  "Fragaria_nilgerrensis",
  "Drymocallis_arguta",
  "Chamaecallis_perpusilloides",
  "Chamaerhodos_erecta",
  "Dasiphora_fruticosa",
  "Potaninia_mongolica",
  "Argentina_anserina",
  "Potentilla_acaulis",
  "Polylepis_tarapacana",
  "Acaena_ovalifolia",
  "Margyricarpus_pinnatus",
  "Sanguisorba_minor",
  "Bencomia_exstipulata",
  "Agrimonia_pilosa_var._pilosa",
  "Rosa_chinensis",
  "Coluria_longifolia",
  "Geum_urbanum",
  "Taihangia_rupestris",
  "Waldsteinia_ternata",
  "Rubus_argutus",
  "Filipendula_ulmaria",
  "Cercocarpus_ledifolius",
  "Purshia_tridentata",
  "Chamaebatia_foliolosa",
  "Dryas_ajanensis_subsp._ajanensis",
  "Neillia_sinensis",
  "Physocarpus_opulifolius",
  "Lyonothamnus_floribundus",
  "Prunus_padus",
  "Aruncus_dioicus",
  "Spiraea_blumei",
  "Prinsepia_uniflora",
  "Exochorda_racemosa_subsp._serratifolia",
  "Oemleria_cerasiformis",
  "Coleogyne_ramosissima",
  "Kerria_japonica",
  "Rhodotypos_scandens",
  "Sorbaria_kirilowii_var._arborea",
  "Gillenia_trifoliata",
  "Lindleya_mespiloides",
  "Vauquelinia_australis",
  "Kageneckia_oblonga",
  "Hesperomeles_cuneata",
  "Phippsiomeles",
  "Dichotomanthes_tristaniicarpa",
  "Cormus_domestica",
  "Chamaemeles_coriacea",
  "Cotoneaster_frigidus",
  "Stranvaesia_nussia",
  "Eriobotrya_japonica",
  "Rhaphiolepis_ferruginea",
  "Sorbus_aucuparia",
  "Pyrus_communis",
  "Karpatiosorbus_bristoliensis",
  "Aria_edulis",
  "Torminalis_glaberrima",
  "Thomsonaria_caloneura",
  "Griffitharia_hemsleyi",
  "Micromeles_alnifolia",
  "Alniaria_alnifolia",
  "Aronia_arbutifolia",
  "Osteomeles_schweriniae",
  "Pseudocydonia_sinensis",
  "Cydonia_oblonga",
  "Chaenomeles_speciosa",
  "Pourthiaea_amphidoxa",
  "Pyracantha_coccinea",
  "Photinia_prunifolia",
  "Malacomeles_denticulata",
  "Amelanchier_laevis",
  "Peraphyllum_ramosissimum",
  "Crataegus_hupehensis",
  "Macromeles_tschonoskii",
  "Malus_domestica",
  "Weniomeles_bodinieri"
)

# 根据不同等级的ILS/IH严苛程度，筛选数据集  --------------------------------
message("4. count triplet dataframe")

## 清理MSC的QT结果
Qtest_clean <- clean_data(Qtest_2)

## 计算物种树的QT矩阵
QT_sp <- get_qt_df(sptree_list,tnames)


## 设置必要参数
alpha_list <- c(0.001, 0.0001, 0.00001, 0.000001)
beta_list <- c(0.01, 0.05)
data_list <- c("sp_data", "red", "yellow", "green")
keep_name_cols <- c("12|34", "13|24", "14|23", "p_T3", "p_star")
keep_sp <- c("12|34", "13|24", "14|23")

## 开始循环
tree_files <- list.files(tree_dir, pattern = "\\.tre$", full.names = TRUE)
message(paste("Found", length(tree_files), "tree files in", tree_dir))


#单个基因树
for(file in tree_files){
  message(paste("Processing:", basename(file)))
  
  for (alpha in alpha_list) {
    for (beta in beta_list) {
      #准备原始的数据集
      blue <- Qtest_clean %>%
        filter(p_T3 >= alpha & p_star <= beta)
      
      red <- Qtest_clean %>%
        filter(p_T3 < alpha & p_star <= beta)
      
      yellow <- Qtest_clean %>%
        filter(p_T3 >= alpha & p_star > beta)
      
      green <- Qtest_clean %>%
        filter(p_T3 < alpha & p_star > beta)
      
      #计算单个基因的QT矩阵
      tree <- read.tree(file)
      tree <- list(tree)
      class(tree) <- "multiPhylo"
      #当前基因树存在的物种个数
	  present_taxa <- taxonNames(tree)
	   # 若少于4个物种则跳过
      if (length(present_taxa) < 4) {
        warning(paste("Tree", basename(file), "has less than 4 taxa, skipped."))
        next
      }

	  # 关键步骤：根据物种树顺序重新排列
      sorted_taxa <- tnames[tnames %in% present_taxa]
	  
      gene_qt <- get_qt_df(tree,sorted_taxa)
   
      #提取单个基因中与物种树拓扑不一致的基因
      name <- tools::file_path_sans_ext(basename(file)) %>% sub("_.*", "",.)
      gene_diff <- gene_qt %>% 
        filter(!gene_qt$comb_all %in% QT_sp$comb_all)
      
      #从这些不一致的基因中，提取属于各个区块的基因
      gene_diff_blue <- gene_diff %>% 
        filter(gene_diff$comb %in% blue$comb)
      
      gene_diff_red <- gene_diff %>% 
        filter(gene_diff$comb %in% red$comb)
      
      gene_diff_yellow <- gene_diff %>% 
        filter(gene_diff$comb %in% yellow$comb)
      
      gene_diff_green <- gene_diff %>% 
        filter(gene_diff$comb %in% green$comb)
      
      #计算当前基因与物种树冲突的quartet数量，以及冲突的quartet中，各个色块的基因数量
      df_num <- data.frame(
        all_quartet=as.numeric(nrow(gene_qt)),
        all_diff=round(as.numeric(nrow(gene_diff)/nrow(gene_qt)),4),
        blue=round(as.numeric(nrow(gene_diff_blue)/nrow(gene_qt)),4),
        red=round(as.numeric(nrow(gene_diff_red)/nrow(gene_qt)),4),
        yellow=round(as.numeric(nrow(gene_diff_yellow)/nrow(gene_qt)),4),
        green=round(as.numeric(nrow(gene_diff_green)/nrow(gene_qt)),4),
        alpha_value=alpha,
        beta_value=beta
      )
      rownames(df_num) <- name
      
      write.csv(df_num, 
                paste0(output_dir, name, "_alpha", alpha, "_beta", beta, "_df_num.csv"), 
                row.names = TRUE)
      
      #将对应的物种对成对提取出来
      # 从 comb 列提取4个物种
      
      all_color_data <- c("gene_diff_red","gene_diff_yellow","gene_diff_green")
      for(col in all_color_data){
        gene_diff_col <- get(col)
        
        # 检查数据框是否为空，如果为空则跳过
        if(nrow(gene_diff_col) <= 1){
          cat("Skipping", col, "- no data available for alpha =", alpha, ", beta =", beta, "\n")
          next  # 直接跳到下一个循环
        }
        
        gene_diff_col <- gene_diff_col %>%
          mutate(
            # 将 comb 列按空格分割成4个物种
            species_list = str_split(comb, "\\s+"),
            taxon1 = map_chr(species_list, ~.x[1]),
            taxon2 = map_chr(species_list, ~.x[2]),
            taxon3 = map_chr(species_list, ~.x[3]),
            taxon4 = map_chr(species_list, ~.x[4])
          )
        
        # 生成物种对并累计计数
        species_pairs <- gene_diff_col %>%
          mutate(row_id = row_number()) %>%
          pmap_dfr(function(taxon1, taxon2, taxon3, taxon4, `12|34`, `13|24`, `14|23`, ...) {
            result <- list()
            
            # 12|34: (taxon1, taxon2) vs (taxon3, taxon4)
            if (!is.na(`12|34`) && `12|34` > 0) {
              result <- c(result, list(
                tibble(
                  pair1 = paste(taxon1, taxon2, sep = "_"),
                  pair2 = paste(taxon3, taxon4, sep = "_"),
                  support = `12|34`,
                  topology = "12|34"
                )
              ))
            }
            
            # 13|24: (taxon1, taxon3) vs (taxon2, taxon4)
            if (!is.na(`13|24`) && `13|24` > 0) {
              result <- c(result, list(
                tibble(
                  pair1 = paste(taxon1, taxon3, sep = "_"),
                  pair2 = paste(taxon2, taxon4, sep = "_"),
                  support = `13|24`,
                  topology = "13|24"
                )
              ))
            }
            
            # 14|23: (taxon1, taxon4) vs (taxon2, taxon3)
            if (!is.na(`14|23`) && `14|23` > 0) {
              result <- c(result, list(
                tibble(
                  pair1 = paste(taxon1, taxon4, sep = "_"),
                  pair2 = paste(taxon2, taxon3, sep = "_"),
                  support = `14|23`,
                  topology = "14|23"
                )
              ))
            }
            
            bind_rows(result)
          })
        
        # 累计计数：对相同的物种对进行求和
        pairs_summary <- species_pairs %>%
          pivot_longer(
            cols = c(pair1, pair2), 
            names_to = "pair_type", 
            values_to = "species_pair"
          ) %>%
          select(species_pair, support, topology) %>%
          mutate(support = as.numeric(support)) %>%
          group_by(species_pair) %>%
          summarise(
            total_support = sum(support), 
            count = n(), 
            .groups = "drop"
          )
        
        # 导出 pairs_summary 为 CSV
        csv_filename <- paste0(output_dir, name, "_", col, "_alpha_", alpha, "_beta_", beta, "_pairs_summary.csv")
        write.csv(pairs_summary, csv_filename, row.names = FALSE)
        
        cat("Exported:", csv_filename, "\n")
        cat("Total species pairs:", nrow(pairs_summary), "\n")
        
        # 6. 创建冲突矩阵
        conflict_matrix <- matrix(
          0,
          nrow = length(species_order),
          ncol = length(species_order),
          dimnames = list(species_order, species_order)
        )
        
        # 7. 填充矩阵
        for (i in seq_len(nrow(pairs_summary))) {
          pair <- pairs_summary$species_pair[i]
          support <- pairs_summary$total_support[i]
          
          sp1 <- NULL
          sp2 <- NULL
          
          for (species_name in species_order) {
            escaped_name <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", species_name)
            if (grepl(paste0("^", escaped_name, "_"), pair)) {
              sp1 <- species_name
              sp2 <- sub(paste0("^", escaped_name, "_"), "", pair)
              if (sp2 %in% species_order) break
            }
          }
          
          if (!is.null(sp1) && !is.null(sp2) && sp1 %in% species_order && sp2 %in% species_order) {
            conflict_matrix[sp1, sp2] <- conflict_matrix[sp1, sp2] + support
            conflict_matrix[sp2, sp1] <- conflict_matrix[sp2, sp1] + support
          }
        }
        
        assign(paste0(name, "_conflict_matrix"), conflict_matrix, envir = .GlobalEnv)
        
        # 导出 conflict_matrix 为 CSV
        
        csv_matrix_filename <- paste0(output_dir, name, "_", col, "_alpha_", alpha, "_beta_", beta, "_conflict_matrix.csv")
        write.csv(conflict_matrix, csv_matrix_filename, row.names = TRUE)
        
        cat("Exported:", csv_matrix_filename, "\n")
        cat("Matrix filled. Non-zero entries:", sum(conflict_matrix > 0), "\n")
        
        # 8. 转为长格式
        matrix_long <- melt(
          conflict_matrix, 
          varnames = c("Species1", "Species2"), 
          value.name = "Support"
        ) %>%
          mutate(Log_Support = log10(Support + 1))
        
        assign(paste0(name, "_matrix_long"), matrix_long, envir = .GlobalEnv)
        
        # 9. 绘制热图
        heatmap_plot <- ggplot(matrix_long, aes(x = Species1, y = Species2, fill = Support)) +
          geom_tile(color = "grey70", linewidth = 0.3) +
          scale_fill_gradientn(
            colors = c("white", "white", "red"),
            values = scales::rescale(c(0, max(matrix_long$Support) * 0.1, max(matrix_long$Support))),
            limits = c(0, max(matrix_long$Support)),
            na.value = "grey95"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 8),
            axis.title = element_blank(),
            panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
          ) +
          coord_fixed() +
          labs(
            title = paste("Species Pair Conflict Matrix -", toupper(name), 
                          "\n(alpha =", alpha, ", beta =", beta, ")"),
            fill = "Support"
          )
        
        # 10. 保存图片
        
        pdf_filename <- paste0(output_dir, name, "_", col, "_alpha_", alpha, "_beta_", beta, "_heatmap_ags353.pdf")
        ggsave(
          filename = pdf_filename, 
          plot = heatmap_plot, 
          width = 14, 
          height = 12
        )
      }
    }
  }
}
