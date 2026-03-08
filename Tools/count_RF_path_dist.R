library(ape)
library(phangorn)
library(ggplot2)

file <- commandArgs(TRUE)

sp <- read.tree(file[1])
gt <- read.tree(file[2])

rf_value <- numeric(length(gt))
path_value <- numeric(length(gt))

for (i in seq_along(gt)) {
  gt1 <- gt[[i]]
  
  common <- intersect(sp$tip.label, gt1$tip.label)
  sp2 <- drop.tip(sp, setdiff(sp$tip.label, common))
  gt2 <- drop.tip(gt1, setdiff(gt1$tip.label, common))
  
  rf_value[i] <- RF.dist(sp2, gt2, check.labels = TRUE,normalize = TRUE)
  path_value[i] <- path.dist(sp2, gt2, check.labels = TRUE, use.weight = FALSE)
}


result <- data.frame(
  RF_Distance = rf_value, 
  Path_Distance = path_value
  )

write.csv(result,"Rf_path_distance_sp_gene.csv",row.names = FALSE,quote = FALSE)

p1 <- ggplot(result, aes(x = Path_Distance, y = RF_Distance)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal() +
  labs(x = "RF Distance", y = "Path Distance",
       title = "Species Tree vs Gene Trees")

############################   MDS   #####################################

# 统一标签：只保留交集部分
fix_labels <- function(sp, gt) {
  common <- intersect(sp$tip.label, gt$tip.label)
  sp_pruned <- keep.tip(sp, common)
  gt_pruned <- keep.tip(gt, common)
  list(sp = sp_pruned, gt = gt_pruned)
}

# 所有树放一起
all_trees <- c(list(sp), gt)

# 计算 RF 距离矩阵
dist_rf <- matrix(0, length(all_trees), length(all_trees))
dist_path <- matrix(0, length(all_trees), length(all_trees))

for (i in 1:(length(all_trees)-1)) {
  for (j in (i+1):length(all_trees)) {
    tmp <- fix_labels(all_trees[[i]], all_trees[[j]])
    dist_rf[i,j] <- RF.dist(tmp$sp, tmp$gt, normalize=TRUE, rooted=FALSE)
    dist_rf[j,i] <- dist_rf[i,j]
    dist_path[i,j] <- path.dist(tmp$sp, tmp$gt, check.labels = TRUE, use.weight = FALSE)
    dist_path[j,i] <- dist_path[i,j]
  }
}
# MDS
RF_ret <- cmdscale(as.dist(dist_rf), k=2)
Path_ret <- cmdscale(as.dist(dist_path), k=2)


#计算RF的均值与标准差，以及MDS矩阵的R方
rf_mean <- mean(dist_rf[1,])
rf_sd <- sd(dist_rf[1,])

dd <- data.frame(mean= rf_mean, sd = rf_sd)
write.csv(dd,"mean_sd_rf.csv",row.names = FALSE, quote = FALSE)


#计算Mean distance to centroid和Mean pairwise distance以及去离群值后的结果
calc_dispersion <- function(coords, remove_outliers = TRUE) {
  # coords: n x k 矩阵，比如 cmdscale 的结果
  # Step1: 点到质心的距离
  centroid <- colMeans(coords)
  dist_to_centroid <- sqrt(rowSums((coords - centroid)^2))
  
  # Step2: 点对距离
  dist_matrix <- dist(coords)

  mdc_norm <- mean(dist_to_centroid)
  mpd_norm <- mean(dist_matrix)
  
  result <- list(
    MDC_norm = mdc_norm,
    MPD_norm = mpd_norm
  )
  
  # Step3: 去掉离群值（如果需要）
  if (remove_outliers) {
    # 使用 median + MAD
    med <- median(dist_to_centroid)
    mad_val <- mad(dist_to_centroid, constant = 1)  # 这里 constant=1 表示原始 MAD
    
    # 定义 Z-score 类似的标准化（基于 median/MAD）
    robust_z <- (dist_to_centroid - med) / mad_val
    keep <- abs(robust_z) <= 3  # 过滤掉 Z > 3 的点
    
    if (sum(keep) > 1) {
      # 重新计算
      centroid2 <- colMeans(coords[keep, , drop = FALSE])
      dist_to_centroid2 <- sqrt(rowSums((coords[keep, ] - centroid2)^2))
      dist_matrix2 <- dist(coords[keep, , drop = FALSE])
      
      result$MDC_norm_no_outlier <- mean(dist_to_centroid2)
      result$MPD_norm_no_outlier <- mean(dist_matrix2)
    } else {
      result$MDC_norm_no_outlier <- NA
      result$MPD_norm_no_outlier <- NA
    }
  }
  
  return(result)
}

res <- calc_dispersion(RF_ret)
res_df <- as.data.frame(t(unlist(res)))

res_path <- calc_dispersion(Path_ret)
res_path_df <- as.data.frame(t(unlist(res_path)))

total_res <- data.frame(
  RF=res_df,
  PATH=res_path_df
)

write.csv(total_res, "dispersion_rf_path.csv", row.names = FALSE, quote = FALSE)

###########################   画图   ######################################
### RF distance
df_rf <- data.frame(
  X = RF_ret[,1],
  Y = RF_ret[,2],
  Type = c("Species tree", rep("Gene tree", length(gt)))
)

text_label <- paste0(
  "MDC = ", round(res_df$MDC_norm_no_outlier, 3), "\n",
  "MPD = ", round(res_df$MPD_norm_no_outlier, 3)
)


p2 <- ggplot(df_rf, aes(x = X, y = Y, color = Type, shape = Type)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Species tree" = "red", "Gene tree" = "blue")) +
  scale_shape_manual(values = c("Species tree" = 19, "Gene tree" = 1)) +
  theme_bw() +
  labs(
    title = "MDS of Species Tree vs Gene Trees (RF Distance)",
    x = "Dimension 1",
    y = "Dimension 2"
  )+
  theme(legend.position = c(0.9,0.9))+
  annotate("text", x = (max(df_rf$X)+min(df_rf$X))/2, y = max(df_rf$Y), 
           label = text_label, hjust = 0, vjust = 1, size = 4)



### PATH distance
df_path <- data.frame(
  X = Path_ret[,1],
  Y = Path_ret[,2],
  Type = c("Species tree", rep("Gene tree", length(gt)))
)

text_label <- paste0(
  "MDC = ", round(res_path$MDC_norm_no_outlier, 3), "\n",
  "MPD = ", round(res_path$MPD_norm_no_outlier, 3)
)


p3 <- ggplot(df_path, aes(x = X, y = Y, color = Type, shape = Type)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Species tree" = "red", "Gene tree" = "blue")) +
  scale_shape_manual(values = c("Species tree" = 19, "Gene tree" = 1)) +
  theme_bw() +
  labs(
    title = "MDS of Species Tree vs Gene Trees (RF Distance)",
    x = "Dimension 1",
    y = "Dimension 2"
  )+
  theme(legend.position = c(0.9,0.9))+
  annotate("text", x = (max(df_path$X)+min(df_path$X))/2, y = max(df_path$Y), 
           label = text_label, hjust = 0, vjust = 1, size = 4)


ggsave("RF_path_rela.pdf", p1, width = 10, height = 8)
ggsave("RF_MDS.pdf", p2, width = 10, height = 8)
ggsave("PATH_MDS.pdf", p3, width = 10, height = 8)