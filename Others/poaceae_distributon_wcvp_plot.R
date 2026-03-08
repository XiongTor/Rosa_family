# !/usr/bin/bash
# Author:Tao Xiong
# Date:2025.05.07
# Description:This script is used to plot the distribution of Lolium perenne and Lolium allgenus.

library(rWCVP)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

data <- read.csv("Lolium_perenne_distribution.csv")
data$latitude <- as.numeric(data$latitude)
data$longitude <- as.numeric(data$longitude)

data2 <- read.csv("Lolium_allgenus.csv")
data2$latitude <- as.numeric(data2$latitude)
data2$longitude <- as.numeric(data2$longitude)


#get the world map
world_map <- map_data("world")

#get the WCVP distribution data
distribution = wcvp_distribution("Lolium",
                                 taxon_rank="species",
                                 extinct=FALSE,
                                 location_doubtful = FALSE)


#清理数据，去除NA
data2_clean <- data2[!is.na(data2$longitude) & !is.na(data2$latitude), ]

# 转化为sf格式
points_sf <- st_as_sf(data2_clean, coords = c("longitude", "latitude"), crs = st_crs(distribution))

# 先把 distribution 转换到投影坐标系（米），适合做距离计算
distribution_proj <- st_transform(distribution, crs = 3857)  # 3857 是 Web Mercator，单位是米
points_proj <- st_transform(points_sf, crs = 3857)
# 扩展区域10km
distribution_buffered <- st_buffer(distribution_proj, dist = 10000)

# 判断点是否在扩展后的范围内
inside <- st_intersects(points_proj, distribution_buffered, sparse = FALSE)

# 找出在区域内外的点
inside_points <- data2_clean[apply(inside, 1, any), ]
outside_points <- data2_clean[!apply(inside, 1, any), ]

#将点映射到不同的wcvp的分区中，即引种还是原生种
inside_proj <- st_transform(
  st_as_sf(inside_points, coords = c("longitude", "latitude"), crs = st_crs(distribution)),
  crs = 3857
)

joined_points <- st_join(inside_proj, distribution_buffered["occurrence_type"])



#画图
col1=c("#995499","#72994c")
# map the WCVP sp dsitribution area to worldmap
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "NA", color = "gray70") +  
  geom_sf(data = distribution, aes(fill = occurrence_type),alpha=0.1,color = "gray60") +  # 分布图层
  geom_sf(data = joined_points, aes(color = occurrence_type), size = 2, alpha = 0.7) +
  scale_fill_manual(values = col1,
                    labels = c("native" = "Native", "introduced" = "Introduced"))+
  scale_color_manual(values = col1) +
  guides(fill = "none")+
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("all_lolium")+
  theme(legend.position = c(0.10, 0.33),
        legend.title = element_blank(),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))


