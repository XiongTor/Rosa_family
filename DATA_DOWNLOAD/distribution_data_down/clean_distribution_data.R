# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025-02-26
# 本脚本用于清理分布数据,但仍有部分要求，要求需要进行清理的数据集必须包含Species、Longitude、Latitude三个字段作为列名

# 加载包
library(dplyr)
library(CoordinateCleaner)

#读取文件
file<-commandArgs(TRUE)
data<-read.csv(file)
data <- as.data.frame(data)


## 保留经纬度非空值
occs.coord <- data[!is.na(data$Latitude) & trimws(as.character(data$Latitude)) != "" &
  !is.na(data$Longitude) & trimws(as.character(data$Longitude)) != "", ]

#CooordinateCleaner包进行数据清理

result <- tryCatch({
  data_problem <- clean_coordinates(
    x = occs.coord,
    lon = "Longitude",
    lat = "Latitude",
    species = "Species",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "seas", "zeros")
  )
  # 如果没有错误，返回结果
  distribution_all_clean <- data[which(data_problem$.summary== "TRUE"),]

  ### 去除重复项
  distribution_duplicate_removal <- distinct(distribution_all_clean,Species,Longitude,Latitude,keep_all = FALSE)
  write.csv(distribution_duplicate_removal,"distribution_duplicate_removal.csv",quote = F,row.names = F)
  return(distribution_all_clean)
}, error = function(e) {
  # 捕获错误信息并存入向量
  error_msg <- conditionMessage(e)
  error_vector <<- c(error_vector, error_msg) # 使用 <<- 修改外部变量
  print(paste("捕获错误：", error_msg))
  return(NULL) # 返回NULL表示出错
})


tt <- capture.output(cat(error_vector))%>%as.data.frame()

tt<-tt[-1,]%>%as.data.frame()

occs.coord<-occs.coord[-tt%tt,]%>%as.data.frame()

occs.coord$Longitude <- as.numeric(occs.coord$Longitude)
occs.coord$Latitude <- as.numeric(occs.coord$Latitude)

data_problem <- clean_coordinates(
    x = occs.coord,
    lon = "Longitude",
    lat = "Latitude",
    species = "Species",
    tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers", "seas", "zeros")
  )
  # 如果没有错误，返回结果
  distribution_all_clean <- data[which(data_problem$.summary== "TRUE"),]

  ### 去除重复项
  distribution_duplicate_removal <- distinct(distribution_all_clean,Species,Longitude,Latitude,keep_all = FALSE)
  write.csv(distribution_duplicate_removal,"distribution_duplicate_removal_out_error.csv",quote = F,row.names = F)
