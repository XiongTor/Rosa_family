# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.03.15
# Description: This script is used to map the Rosa subfamilies coordinate information to the world map.

library(ggplot2)
library(maps)
library(dplyr)

data <- read.csv("Final_all_distribution_cleaned.csv")
data2 <- read.csv("native/All_Native.csv")

rosa <- read.csv("rosa_species_subfam_inf.csv")

Rosoideae <- rosa %>% 
  filter(Subfamilies=="Rosoideae")

Dryadoideae <- rosa %>% 
  filter(Subfamilies=="Dryadoideae")

Amygdaloideae <- rosa %>% 
  filter(Subfamilies=="Amygdaloideae")

Rosoideae_distribution <- data2 %>% 
  filter(species_name%in%Rosoideae$Species)

Dryadoideae_distribution <- data2 %>% 
  filter(species_name%in%Dryadoideae$Species)

Amygdaloideae_distribution <- data2 %>% 
  filter(species_name%in%Amygdaloideae$Species)


# 绘制世界地图
world_map <- map_data("world")
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "lightblue", color = "white") +
  geom_point(data = Rosoideae_distribution, aes(x = combined_lon, y = combined_lat), color = "red", size = 0.5) +
  theme_minimal() +
  labs(title = "Rosoideae_native")
