# 地理分布数据下载
library(dplyr)
library(sf)
library(spData)
library(tmap)
library(purrr)
library(readr)
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid
library(maptools)
library(raster)
library(ggplot2)
library(ggcorrplot)
library(CoordinateCleaner)
library(zoo)
library(ridigbio)
library(readxl)



## 下载gbif数据（提前准备物种名录）
file_url <- read.csv("rosa_species_subfam_inf.csv", header=TRUE, sep=";")

#设置请求的服务器IP，即VPN挂的地址
Sys.setenv(ALL_PROXY = "127.0.0.1:7890")

#下载GBIF数据
library(rgbif)

for (i in file_url$Species) {
  tryCatch({
    print(paste("Downloading GBIF data for species:", i))
    
    # 1. 使用 occ_search 并移除已废弃的 return = "data" 参数
    item_gbif <- occ_search(
      scientificName = i,
      limit = 90000,
      hasCoordinate = TRUE
    )
    
    # 2. 检查数据是否存在
    if (!is.null(item_gbif$data) && nrow(item_gbif$data) > 0) {
      # 3. 提取实际的 data 部分
      item_gbif_data <- item_gbif$data
      
      # 4. 检查需要的列是否存在
      required_cols <- c("species", "decimalLongitude", "decimalLatitude")
      if (all(required_cols %in% colnames(item_gbif_data))) {
        # 5. 筛选有效坐标（排除 0,0 坐标）
        item_gbif_reduced_cleaned <- item_gbif_data[
          item_gbif_data$decimalLatitude != 0 & 
          item_gbif_data$decimalLongitude != 0,
          required_cols
        ]
        
        # 6. 创建结果数据框（更简洁的方式）
        combined <- data.frame(
          species_name = item_gbif_reduced_cleaned$species,
          lon = item_gbif_reduced_cleaned$decimalLongitude,
          lat = item_gbif_reduced_cleaned$decimalLatitude
        )
        
        # 7. 这里可以添加保存或合并结果的代码
        # 比如写入文件或累加到总数据框
        write.table(combined, file = paste0("./results/", gsub(" ", "_", i), ".csv", sep = ""), sep = ",", quote = FALSE, row.names = FALSE)
        
        print(paste("Successfully processed", nrow(combined), "records for", i))
      } else {
        print(paste("Missing required columns for", i))
      }
    } else {
      print(paste("No data found for", i))
    }
  }, error = function(e) {
    message(paste("Error processing", i, ":", conditionMessage(e)))
  })
}



#===============================================================================


## 下载idigbio数据
library("ridigbio")
species_list <- read.csv("rosa_species_subfam_inf.csv", header=TRUE, sep=";")
species_names <- species_list$Species
all_results <- list()  # 初始化存储结果的列表 

already_down <- data.frame(genus=names(all_results))

remaining_sp <- species_names[!species_list$Species %in% already_down$genus]


for (species in remaining_sp) {
  rq <- list(scientificname = species)  # 创建查询条件
  results <- idig_search_records(rq)    # 获取结果
  all_results[[species]] <- results     # 将结果存储在列表中
}


combined_results <- do.call(rbind, all_results)  
#过滤一下经纬度为空的数值
final_idig <- combined_results %>% 
  filter(combined_results$geopoint.lon!="NA" &combined_results$geopoint.lat!="NA")

#清理数据
data_problem <- clean_coordinates(x = final_idig,   
                                  lon = "geopoint.lon", 
                                  lat = "geopoint.lat",
                                  species = "scientificname",
                                  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers",
                                            "seas", "zeros")) 


summary(data_problem)
data_problematic <- data[which(data_problem$.summary== "FALSE"),]
distribution_all_clean <- data[which(data_problem$.summary== "TRUE"),]

write.csv(final_idig,file="idigbio_rosaceae_2025.2.23.csv")
