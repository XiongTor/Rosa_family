library(ape)
library(phangorn)
library(MSCquartets)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(doParallel)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(viridis)
library(reshape2)

# 1. 读取数据并计算物种树种的四连体分布情况
# orthofinder  
# load("12-MSCquartet/orthofinder/quartet_analysis_results.RData")
# data <- as.data.frame(Qtest_2)

## ags353
load("12-MSCquartet/ags353/quartet_analysis_results.RData")
data <- as.data.frame(Qtest_2)
target_cols <- c("13|24", "12|34", "14|23")
data[, target_cols] <- data[, target_cols] / rowSums(data[, target_cols])

# 查看结果前几行以确认
head(data[, target_cols])

rm(Qtest_2,pTable,pTable_2)
## sptree
sptree=read.tree("12-MSCquartet/ags353/rosa_ags353_treeshrink_sp_rt_oneoutg_final.tre")
sptree2 <- write.tree(sptree)
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"
### 得到数据集中所有物种名字
tnames=taxonNames(sptree_list)
## 统计基因树中所有可能的quartet的三种拓扑出现次数
QT_sptree <- quartetTable(sptree_list,tnames)
data_sptree <- as.data.frame(QT_sptree)
rm(QT_sptree)


# 2. 计算简化相关数据
data_list <- c("sp_data","red","yellow","green")
keep_name_cols <- c("12|34", "13|24", "14|23", "p_T3", "p_star")
keep_sp <- c("12|34", "13|24", "14|23")

alpha_list <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001)
beta_list <- c(0.01, 0.05)
for(alpha in alpha_list){
  for(beta in beta_list){
    alpha <- 0.01
    beta <- 0.05
    
    sp_data <- data_sptree %>% 
      mutate(across(-all_of(keep_sp), ~ifelse(. == 1, cur_column(), NA)))
    
    # blue <- data %>%
    #   filter(data$p_T3 >= alpha & data$p_star<= beta) %>% 
    #   mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
    
    red <- data %>%
      filter(p_T3 < alpha & p_star <= beta) %>% 
      mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
    
    yellow <- data %>%
      filter(p_T3 >= alpha & p_star > beta) %>%
      mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
    
    green <- data %>%
      filter(p_T3 < alpha & p_star > beta) %>% 
      mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
    
    
    for(name in data_list){
      # 获取当前数据框
      current_data <- get(name)
      
      # 处理数据
      sorted_data <- current_data %>%
        {
          t(apply(., 1, function(row) {
            non_na <- row[!is.na(row)]
            c(non_na, rep(NA, length(row) - length(non_na)))
          }))
        } %>%
        as.data.frame() %>% 
        select(where(~!all(is.na(.))))
      
      # 重命名列
      # 根据 name 选择列名向量
      if (name == "sp_data") {
        name_cols <- keep_sp
      } else {
        name_cols <- keep_name_cols
      }
      
      new_colnames <- c(
        paste0("taxon", 1:4),
        name_cols
      )
      
      colnames(sorted_data) <- new_colnames
      
      # 保存回环境中（新变量名）
      assign(paste0(name, "_sorted"), sorted_data, envir = .GlobalEnv)
    }
    
    
    # 3. 保留各个色块数据框中，与物种树拓扑不一致的四连体的数值，若一致则记为0视为不参与贡献-----------------------------------------------
    ## 合并数据集合，将不同列的名字合并到一起作为一个字符
    # 定义数据集列表
    data_list2 <- c("red", "yellow")
    
    # 为物种树数据添加查找键（在循环外，只需做一次）
    sp_data_with_key <- sp_data_sorted %>%
      mutate(
        quartet_key = paste(taxon1, taxon2, taxon3, taxon4, sep = "|"),
        sp_topology = case_when(
          `12|34` == 1 ~ "12|34",
          `13|24` == 1 ~ "13|24",
          `14|23` == 1 ~ "14|23"
        )
      ) %>%
      select(quartet_key, sp_topology)
    
    # 主循环开始---------------------------------------------------
    for(name in data_list2) {
      
      cat("\n========== Processing:", name, "==========\n")
      
      # 1. 获取当前数据集
      current_data <- get(paste0(name, "_sorted"))
      
      # 2. 添加quartet_key
      data_with_key <- current_data %>%
        mutate(quartet_key = paste(taxon1, taxon2, taxon3, taxon4, sep = "|"))
      
      # 保存：name_with_key
      assign(paste0(name, "_with_key"), data_with_key, envir = .GlobalEnv)
      
      # 3. 与物种树比对，保留不一致的拓扑
      data_filtered <- data_with_key %>%
        left_join(sp_data_with_key, by = "quartet_key") %>%
        mutate(
          # 统一转换拓扑列为数值型
          across(c(`12|34`, `13|24`, `14|23`), as.numeric),
          # 过滤掉与物种树一致的拓扑
          `12|34` = ifelse(!is.na(sp_topology) & sp_topology == "12|34", 0, `12|34`),
          `13|24` = ifelse(!is.na(sp_topology) & sp_topology == "13|24", 0, `13|24`),
          `14|23` = ifelse(!is.na(sp_topology) & sp_topology == "14|23", 0, `14|23`)
        ) %>%
        select(-quartet_key, -sp_topology)
      
      # 保存：name_filtered
      assign(paste0(name, "_filtered"), data_filtered, envir = .GlobalEnv)
      
      # 4. 将四元组转化为物种对
      species_pairs <- data_filtered %>%
        mutate(row_id = row_number()) %>%
        pmap_dfr(function(taxon1, taxon2, taxon3, taxon4, `12|34`, `13|24`, `14|23`, ...) {
          result <- list()
          
          # 12|34 拓扑
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
          
          # 13|24 拓扑
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
          
          # 14|23 拓扑
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
      
      # 保存：name_species_pairs
      assign(paste0(name, "_species_pairs"), species_pairs, envir = .GlobalEnv)
      
      # 5. 合并物种对并统计
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
      
      # 保存：name_pairs_summary
      assign(paste0(name, "_pairs_summary"), pairs_summary, envir = .GlobalEnv)
      
      cat("Total species pairs:", nrow(pairs_summary), "\n")
      
      # 6. 创建物种对矩阵
      # First, define your desired order of species
      tnames <- c(
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
      
      # Create the conflict matrix with the specified order
      conflict_matrix <- matrix(
        0,
        nrow = length(tnames),
        ncol = length(tnames),
        dimnames = list(tnames, tnames)
      )
      
      
      # 7. 填充矩阵
      for (i in seq_len(nrow(pairs_summary))) {
        pair <- pairs_summary$species_pair[i]
        support <- pairs_summary$total_support[i]
        
        sp1 <- NULL
        sp2 <- NULL
        
        # 智能分割物种对
        for (species_name in tnames) {
          escaped_name <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", species_name)
          if (grepl(paste0("^", escaped_name, "_"), pair)) {
            sp1 <- species_name
            sp2 <- sub(paste0("^", escaped_name, "_"), "", pair)
            if (sp2 %in% tnames) break
          }
        }
        
        # 填充矩阵（对称）
        if (!is.null(sp1) && !is.null(sp2) && sp1 %in% tnames && sp2 %in% tnames) {
          conflict_matrix[sp1, sp2] <- conflict_matrix[sp1, sp2] + support
          conflict_matrix[sp2, sp1] <- conflict_matrix[sp2, sp1] + support
        }
      }
      
      # 保存：name_conflict_matrix
      assign(paste0(name, "_conflict_matrix"), conflict_matrix, envir = .GlobalEnv)
      
      cat("Matrix filled. Non-zero entries:", sum(conflict_matrix > 0), "\n")
      
      # 8. 转为长格式
      matrix_long <- melt(
        conflict_matrix, 
        varnames = c("Species1", "Species2"), 
        value.name = "Support"
      ) %>%
        mutate(
          Log_Support = log10(Support + 1)  # +1 避免log(0)
        )
      
      # 保存：name_matrix_long
      assign(paste0(name, "_matrix_long"), matrix_long, envir = .GlobalEnv)
      
      # 9. 绘制热图
      heatmap_plot <- ggplot(matrix_long, aes(x = Species1, y = Species2, fill = Support)) +
        geom_tile(color = "grey70", linewidth = 0.3) +
        scale_fill_gradientn(
          colors = c("white", "white", "red"),
          values = scales::rescale(c(0, max(matrix_long$Support)*0.1, max(matrix_long$Support))),  # 自定义断点
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
          title = paste("Species Pair Conflict Matrix -", toupper(name)),
          fill = "Support"
        )
      
      # 10. 显示并保存图片
      print(heatmap_plot)
      ggsave(
        filename = paste0(name,"alpha_",alpha,"_beta_",beta, "_heatmap_ags353.pdf"), 
        plot = heatmap_plot, 
        width = 14, 
        height = 12
      )
      
      cat("Saved outputs for", name, "\n")
      cat("==========================================\n\n")
    }
  }
}



