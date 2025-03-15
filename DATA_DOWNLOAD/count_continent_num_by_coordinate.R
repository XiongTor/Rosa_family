# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.03.15
# Description: This script is used to calculate the number of native and introduced species in each continent.

# 加载包
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggplot2)
library(maps)

#####################统计native分布的数据#####################################

data <- read.csv("native/All_Native.csv")

# 转换为sf空间对象（使用WGS84坐标系）
sf <- st_as_sf(data, coords = c("combined_lon", "combined_lat"), crs = 4326)

# 下载大洲多边形数据（scale可选50m或10m，10m更精确但体积更大）
continents <- ne_countries(scale = "medium", returnclass = "sf") %>%
  group_by(continent) %>%
  summarise()  # 合并同一大洲的多边形

# 空间连接：为每个点找到相交的大洲
result <- st_join(sf, continents, join = st_intersects)

# 转换为普通数据框并清理结果
final_df <- result %>%
  st_drop_geometry() %>%
  select(species_name,continent) %>% 
  filter(continent!="NA") %>%
  distinct()



####################统计带引种信息的数据#####################################
#合并数据
data <- read.csv("Final_all_distribution_cleaned.csv")

# 转换为sf空间对象（使用WGS84坐标系）
sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

# 下载大洲多边形数据（scale可选50m或10m，10m更精确但体积更大）
continents <- ne_countries(scale = "medium", returnclass = "sf") %>%
  group_by(continent) %>%
  summarise()  # 合并同一大洲的多边形

# 空间连接：为每个点找到相交的大洲
result <- st_join(sf, continents, join = st_intersects)

# 转换为普通数据框并清理结果
final_df2 <- result %>%
  st_drop_geometry() %>%
  select(taxon_name,continent) %>% 
  filter(continent!="NA") %>%
  distinct()

write.csv(final_df,"earth_number/Final_sp_native_earth.csv",row.names = F)
write.csv(final_df2,"earth_number/Final_sp_native_introduced.csv",row.names = F)

#统计属级别的数据
final_df$genus <-sapply(strsplit(final_df$species_name," "),'[',1)

final_df_genus <- final_df %>% 
  select(genus,continent) %>% 
  distinct() %>% 
  group_by(genus) %>% 
  summarise(across(everything(),~paste(unique(.),collapse = ",")))

write.csv(final_df_genus,"earth_number/Final_genus_native_earth.csv",row.names = F)

#=======================================
final_df2$genus <-sapply(strsplit(final_df2$taxon_name," "),'[',1)

final_df2_genus <- final_df2 %>% 
  select(genus,continent) %>% 
  distinct() %>% 
  group_by(genus) %>% 
  summarise(across(everything(),~paste(unique(.),collapse = ",")))

write.csv(final_df_genus,"earth_number/Final_genus_native_introduced_earth.csv",row.names = F)

  



###################统计原种和引种在分布大洲上的差别###########################
# 定义增强函数：返回差异数、相同项、不同项
calculate_diff_detail <- function(cont1, cont2) {
  # 清洗和标准化
  process <- function(s) {
    unique(trimws(tolower(unlist(strsplit(s, ",")))))
  }
  
  set1 <- process(cont1)
  set2 <- process(cont2)
  
  # 计算交集（相同项）
  common <- intersect(set1, set2)
  # 计算对称差集（不同项）
  diff <- union(setdiff(set1, set2), setdiff(set2, set1))
  
  # 转换为逗号分隔字符串（若空则标记为NA）
  list(
    diff_count = length(diff),
    common_continents = ifelse(length(common) == 0, NA, paste(tools::toTitleCase(common), collapse = ", ")),
    diff_continents = ifelse(length(diff) == 0, NA, paste(tools::toTitleCase(diff), collapse = ", "))
  )
}

# 预处理：合并每个物种的分布字符串
species_cont <- final_df_genus %>%
  group_by(genus) %>%
  summarise(cont_species = paste(continent, collapse = ","))

taxon_cont <- final_df2_genus %>%
  group_by(genus) %>%
  summarise(cont_taxon = paste(continent, collapse = ","))

# 执行比较
result <- inner_join(
  species_cont,
  taxon_cont,
  by = "genus"
) %>% 
  rowwise() %>%
  mutate(
    result = list(calculate_diff_detail(cont_species, cont_taxon))
  ) %>%
  unnest_wider(result) %>%  # 展开列表为多列
  select(
    genus,
    diff_count,
    common_continents,
    diff_continents
  )



# 输出结果
write.csv(result,"earth_number/Final_genus_diff_num.csv",row.names = F)







