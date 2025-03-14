# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.03.14
# Description: This script is used to combine different datasets to be one dataset. And then to change the species name to the accpeted name.
#Then to get the native distribution by WCVP

library(dplyr)
library(rWCVP)
library(rgbif)
library(tidyverse)
library(sf)


#合并数据
cvh_nsii <- read.csv("winscp/Final_cleaned_cvh_nsii.csv")
gbif <- read.csv("winscp/Final_GBIF_cleaned_distribution.csv")
idigbio <- read.csv("winscp/Final_idigbio_cleaned_distribution.csv")
sun <- read.csv("winscp/Final_sun_rosaceae_distribution.csv")


cvh_nsii <- cvh_nsii %>% 
  select(Species,Longitude,Latitude)

cvh_nsii$type <- "cvh_nsii"

gbif <- gbif %>% 
  select(-keep_all)
gbif$type <- "gbif"

idigbio <- idigbio %>% 
  select(-keep_all)
idigbio$type <- "idigbio"

sun <- sun %>% 
  select(-keep_all)
sun$type <- "sun"

data <- rbind(cvh_nsii,gbif,idigbio,sun)

#调整数据格式
data <- data %>% 
  mutate(Species=str_replace(Species,"^[a-z]",function(x) toupper(x)))

#获取物种名单

sp_list <- data$Species %>% sort() %>% unique() %>% as.data.frame()
colnames(sp_list) <- "Species"


######################################## 名称标准化 ##############################################

globalWoodiness_matches <- wcvp_match_names(
  sp_list,
  name_col = "Species",
  fuzzy = FALSE,
  progress_bar = FALSE
)

manually_resolved <- filter(
  globalWoodiness_matches,
  wcvp_status != "Illegitimate" &
    wcvp_status != "Misapplied" &
    wcvp_status != "Invalid" &
    wcvp_status != "NA"
)


rWCVP <- rWCVPdata::wcvp_names
rWCVP$plant_name_id <- as.double(rWCVP$plant_name_id)

accepted_matches <- manually_resolved %>%
  left_join(rWCVP, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  mutate(
    keep = case_when(
      taxon_status == "Accepted" & (wcvp_status != "Synonym" | wcvp_homotypic) ~
        "Matched to an accepted name",
      TRUE ~ "Not matched to an accepted name"
    )
  )


globalWoodiness_accepted <- accepted_matches %>%
  filter(accepted_plant_name_id %in% accepted_matches$wcvp_accepted_id) %>%
  group_by(Species) %>%
  filter(
    sum(wcvp_status == "Unplaced") == n() |
      any(wcvp_status != "Unplaced")
  ) %>%
  ungroup() %>%
  select(Species, taxon_name, wcvp_status) %>%
  filter(taxon_name!="NA") %>% 
  distinct(Species,.keep_all=T)


#根据标准化后的物种名，选出能被wcvp匹配上的物种
data2 <- data %>% 
  filter(Species%in%globalWoodiness_accepted$Species)




#将原来的名字替换为标准名称

data2 <- merge(data2,globalWoodiness_accepted,by="Species",all.x = T)

write.csv(data2,"Final_all_distribution_cleaned.csv",row.names = F)


#################################根据wcvp结合kew的自然种和引种信息，筛选出自然种######################

#获得物种列表
sp_list <- data2$taxon_name %>% sort() %>% unique()

#创建空数据框
all_filtered_df <- data.frame()

for(name in sp_list){
  
#提取出需要的数据集
  data_tmp <- data2 %>% 
    filter(taxon_name == name)
  
#将数据集转换为空间数据集
  occs <- data_tmp %>% 
    st_as_sf(coords=c("Longitude", "Latitude"), crs=st_crs(4326))
  
#下载物种的原生分布区范围,并判断是否能捕获到原产地信息
  native_range <- tryCatch({
    wcvp_distribution(name,
                      taxon_rank = "species",
                      introduced = FALSE,
                      extinct = FALSE,
                      location_doubtful = FALSE)
  }, error = function(e) {
    # Append species name and reason to file if no distribution data is available
    message(paste("### Skipping species:", name, "- Reason: No distribution data or name mismatch ###"))
    write(paste("### Skipping species:", name, "- Reason: No distribution data or name mismatch ###"), file = "./missing_native_ranges.txt", append = TRUE, sep = "\n")
    return(NULL)
  })
  
#如果没有找到原产地信息则直接结束循环
  if (is.null(native_range) || nrow(native_range) == 0) next
  
#限制原产地分布范围1km外边界的点都可接收
  occs$native <- st_intersects(occs, st_union(native_range), sparse = FALSE)[,1]
  buffered_dist <- native_range %>%
    st_union() %>%
    st_buffer(0.009)
  occs$native_buffer <- st_intersects(occs, buffered_dist, sparse = FALSE)[,1]
  
#过滤数据
  filtered <- occs %>% filter(native_buffer)
  
#检查是否存在空值的现象
  if (nrow(filtered) == 0) {
    message(paste("### Species", name, "is missing geometry information ###"))
    next
  }
  
#确认filter是否为sf数据格式
  if (!inherits(filtered, "sf") || is.null(st_geometry(filtered))) {
    message(paste("### Species", name, "is missing geometry information ###"))
    next
  }
  
#提取经纬度
  coords <- st_coordinates(filtered)
  filtered_df <- data.frame(
    species_name = filtered$taxon_name,
    combined_lon = coords[, "X"],
    combined_lat = coords[, "Y"]
  )
  
  write.csv(filtered_df, paste0("./native/", gsub(" ", "_", name), "_native.csv"), row.names = FALSE)
  
  # 合并数据
  all_filtered_df <- bind_rows(all_filtered_df, filtered_df)
}

write.csv(all_filtered_df, "./native/All_Native.csv", row.names = FALSE)

message("### done,missing dats save missing_native_ranges.txt。 ###")










