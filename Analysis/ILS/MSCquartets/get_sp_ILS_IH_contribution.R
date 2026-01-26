library(ape)
library(phangorn)
library(MSCquartets)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(viridis)
library(reshape2)

# 1. 读取数据并计算物种树种的四连体分布情况
load("quartet_analysis_results.RData")
data <- as.data.frame(Qtest_2)
target_cols <- c("13|24", "12|34", "14|23")
data[, target_cols] <- data[, target_cols] / rowSums(data[, target_cols])

# 查看结果前几行以确认
head(data[, target_cols])

rm(Qtest_2, pTable, pTable_2)

## sptree
sptree <- read.tree("rosa_orthofinder_MO_treeshrink_sp_rt_oneoutg_final.tre")
sptree2 <- write.tree(sptree)
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"

### 得到数据集中所有物种名字
tnames <- taxonNames(sptree_list)

## 统计基因树中所有可能的quartet的三种拓扑出现次数
QT_sptree <- quartetTable(sptree_list, tnames)
data_sptree <- as.data.frame(QT_sptree)
rm(QT_sptree)

# 2. 设置参数
data_list <- c("sp_data", "red", "yellow", "green")
keep_name_cols <- c("12|34", "13|24", "14|23", "p_T3", "p_star")
keep_sp <- c("12|34", "13|24", "14|23")

alpha_list <- c(0.001, 0.0001, 0.00001, 0.000001)
beta_list <- c(0.01, 0.05)

# 定义物种顺序
species_order <- c(
  "Weniomeles_bodinieri", "Stranvaesia_nussia", "Cormus_domestica",
  "Macromeles_tschonoskii", "Malus_domestica", "Pyrus_communis",
  "Micromeles_alnifolia", "Alniaria_alnifolia", "Aria_edulis",
  "Karpatiosorbus_bristoliensis", "Griffitharia_hemsleyi", "Thomsonaria_caloneura",
  "Torminalis_glaberrima", "Aronia_arbutifolia", "Pseudocydonia_sinensis",
  "Cydonia_oblonga", "Chaenomeles_speciosa", "Pourthiaea_amphidoxa",
  "Osteomeles_schweriniae", "Pyracantha_coccinea", "Eriobotrya_japonica",
  "Rhaphiolepis_ferruginea", "Sorbus_aucuparia", "Cotoneaster_frigidus",
  "Photinia_prunifolia", "Dichotomanthes_tristaniicarpa", "Chamaemeles_coriacea",
  "Amelanchier_laevis", "Malacomeles_denticulata", "Peraphyllum_ramosissimum",
  "Phippsiomeles", "Crataegus_hupehensis", "Hesperomeles_cuneata",
  "Vauquelinia_australis", "Kageneckia_oblonga", "Lindleya_mespiloides",
  "Gillenia_trifoliata", "Prinsepia_uniflora", "Exochorda_racemosa_subsp_serratifolia",
  "Oemleria_cerasiformis", "Coleogyne_ramosissima", "Kerria_japonica",
  "Rhodotypos_scandens", "Sorbaria_kirilowii_var_arborea", "Lyonothamnus_floribundus",
  "Prunus_padus", "Aruncus_dioicus", "Spiraea_blumei",
  "Neillia_sinensis", "Physocarpus_opulifolius", "Purshia_tridentata",
  "Chamaebatia_foliolosa", "Cercocarpus_ledifolius", "Dryas_ajanensis_subsp_ajanensis",
  "Drymocallis_arguta", "Chamaecallis_perpusilloides", "Chamaerhodos_erecta",
  "Dasiphora_fruticosa", "Potaninia_mongolica", "Sibbaldianthe_bifurca",
  "Sibbaldia_parviflora", "Alchemilla_faeroensis", "Comarum_palustre",
  "Fragaria_nilgerrensis", "Argentina_anserina", "Potentilla_acaulis",
  "Acaena_ovalifolia", "Margyricarpus_pinnatus", "Polylepis_tarapacana",
  "Bencomia_exstipulata", "Sanguisorba_minor", "Agrimonia_pilosa_var_pilosa",
  "Rosa_chinensis", "Taihangia_rupestris", "Waldsteinia_ternata",
  "Coluria_longifolia", "Geum_urbanum", "Rubus_argutus",
  "Filipendula_ulmaria", "Outgroup"
)

# 创建输出目录
if (!dir.exists("output_csvs")) {
  dir.create("output_csvs")
}

# 主循环
for (alpha in alpha_list) {
  for (beta in beta_list) {
    
    cat("\n\n########## Processing alpha =", alpha, ", beta =", beta, "##########\n")
    
    # 添加错误处理
    tryCatch({
      
      # 准备数据
      sp_data <- data_sptree %>% 
        mutate(across(-all_of(keep_sp), ~ifelse(. == 1, cur_column(), NA)))
      
      red <- data %>%
        filter(p_T3 < alpha & p_star <= beta) %>% 
        mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
      
      yellow <- data %>%
        filter(p_T3 >= alpha & p_star > beta) %>%
        mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
      
      green <- data %>%
        filter(p_T3 < alpha & p_star > beta) %>% 
        mutate(across(-all_of(keep_name_cols), ~ifelse(. == 1, cur_column(), NA)))
      
      # 排序和重命名
      for (name in data_list) {
        current_data <- get(name)
        
        # 检查数据是否为空，如果为空就跳过
        if (nrow(current_data) == 0) {
          cat("  ⊘ Skipping", name, "(no data)\n")
          next
        }
        
        sorted_data <- current_data %>%
          {
            t(apply(., 1, function(row) {
              non_na <- row[!is.na(row)]
              c(non_na, rep(NA, length(row) - length(non_na)))
            }))
          } %>%
          as.data.frame() %>% 
          select(where(~!all(is.na(.))))
        
        # 检查排序后的列数，如果为0就跳过
        if (ncol(sorted_data) == 0) {
          cat("  ⊘ Skipping", name, "(no valid columns after processing)\n")
          next
        }
        
        if (name == "sp_data") {
          name_cols <- keep_sp
        } else {
          name_cols <- keep_name_cols
        }
        
        new_colnames <- c(paste0("taxon", 1:4), name_cols)
        
        # 确保列数匹配，如果不匹配就跳过
        if (ncol(sorted_data) != length(new_colnames)) {
          cat("  ⊘ Skipping", name, "(column mismatch: expected", length(new_colnames), ", got", ncol(sorted_data), ")\n")
          next
        }
        
        colnames(sorted_data) <- new_colnames
        assign(paste0(name, "_sorted"), sorted_data, envir = .GlobalEnv)
        cat("  ✓ Processed", name, "\n")
      }
      
      # 准备物种树比对数据
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
      
      # 处理每个数据集
      data_list2 <- c("red", "yellow")
      
      for (name in data_list2) {
        
        cat("\n========== Processing:", name, "==========\n")
        
        # 检查这个数据集是否存在且已处理
        sorted_name <- paste0(name, "_sorted")
        if (!exists(sorted_name)) {
          cat("  ⊘ Skipping", name, "(not processed in previous step)\n")
          cat("==========================================\n\n")
          next
        }
        
        # 1. 获取当前数据集
        current_data <- get(sorted_name)
        
        # 检查是否有数据
        if (nrow(current_data) == 0) {
          cat("  ⊘ Skipping", name, "(no data available)\n")
          cat("==========================================\n\n")
          next
        }
        
        # 2. 添加quartet_key
        data_with_key <- current_data %>%
          mutate(quartet_key = paste(taxon1, taxon2, taxon3, taxon4, sep = "|"))
        
        assign(paste0(name, "_with_key"), data_with_key, envir = .GlobalEnv)
        
        # 3. 与物种树比对
        data_filtered <- data_with_key %>%
          left_join(sp_data_with_key, by = "quartet_key") %>%
          mutate(
            across(c(`12|34`, `13|24`, `14|23`), as.numeric),
            `12|34` = ifelse(!is.na(sp_topology) & sp_topology == "12|34", 0, `12|34`),
            `13|24` = ifelse(!is.na(sp_topology) & sp_topology == "13|24", 0, `13|24`),
            `14|23` = ifelse(!is.na(sp_topology) & sp_topology == "14|23", 0, `14|23`)
          ) %>%
          select(-quartet_key, -sp_topology)
        
        assign(paste0(name, "_filtered"), data_filtered, envir = .GlobalEnv)
        
        # 4. 转化为物种对
        species_pairs <- data_filtered %>%
          mutate(row_id = row_number()) %>%
          pmap_dfr(function(taxon1, taxon2, taxon3, taxon4, `12|34`, `13|24`, `14|23`, ...) {
            result <- list()
            
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
        
        assign(paste0(name, "_species_pairs"), species_pairs, envir = .GlobalEnv)
        
        # 5. 统计物种对
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
        
        assign(paste0(name, "_pairs_summary"), pairs_summary, envir = .GlobalEnv)
        
        # 导出 pairs_summary 为 CSV
        csv_filename <- paste0("output_csvs/", name, "_alpha_", alpha, "_beta_", beta, "_pairs_summary.csv")
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
        csv_matrix_filename <- paste0("output_csvs/", name, "_alpha_", alpha, "_beta_", beta, "_conflict_matrix.csv")
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
        pdf_filename <- paste0(name, "_alpha_", alpha, "_beta_", beta, "_heatmap_ags353.pdf")
        ggsave(
          filename = pdf_filename, 
          plot = heatmap_plot, 
          width = 14, 
          height = 12
        )
        
        cat("Saved outputs for", name, "\n")
        cat("==========================================\n\n")
      }
      
    }, error = function(e) {
      cat("ERROR in alpha =", alpha, ", beta =", beta, ":\n")
      cat(conditionMessage(e), "\n")
      cat("Continuing to next iteration...\n\n")
    })
  }
}

cat("\n\n########## All processing completed! ##########\n")
cat("CSV files saved in 'output_csvs/' directory\n")
cat("PDF files saved in current directory\n")