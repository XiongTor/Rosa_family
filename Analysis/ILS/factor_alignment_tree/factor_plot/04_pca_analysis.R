# !/usr/bin/Rscript
# Author: XiongTao
# Description: PCA analysis

library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)

# Load data
merged_orthofinder <- readRDS("merged_orthofinder.rds")
merged_chloroplast <- readRDS("merged_chloroplast.rds")

# Add type labels
merged_orthofinder$type <- "nuclear"
merged_chloroplast$type <- "chloroplast"

# Combine datasets
df_all <- bind_rows(merged_orthofinder, merged_chloroplast)

# Select properties for PCA
props <- c("GC_content", "Proportion_variable_sites", "Proportion_parsimony_informative",
           "Mean_evolutionary_rate", "RCV")

# Perform PCA
X <- scale(df_all[, props])


pca <- prcomp(X, center = TRUE, scale. = TRUE)

saveRDS(df_all, file = "df_all.rds")
saveRDS(pca, file = "pca.rds")

# Extract PCA results
pca_df <- data.frame(
  Alignment_name = df_all$Alignment_name,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  type = df_all$type
)

# Merge PC scores back to original data
pc_map <- pca_df[, c("Alignment_name", "PC1", "PC2", "PC3")]

merged_orthofinder <- merge(merged_orthofinder, pc_map, by = "Alignment_name", all.x = TRUE)
merged_chloroplast <- merge(merged_chloroplast, pc_map, by = "Alignment_name", all.x = TRUE)

saveRDS(merged_orthofinder, "merged_orthofinder.rds")
saveRDS(merged_chloroplast, "merged_chloroplast.rds")

# 2D PCA plots
p1 <- ggplot(pca_df, aes(PC1, PC2, color = type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.15, 0.9), legend.background = element_blank())

p2 <- ggplot(pca_df, aes(PC1, PC3, color = type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.15, 0.9), legend.background = element_blank())

p3 <- ggplot(pca_df, aes(PC2, PC3, color = type)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 1) +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.15, 0.9), legend.background = element_blank())

pdf("pca.pdf", width = 18, height = 6)
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

# 3D PCA plot
pca_df_3d <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  type = df_all$type
)

plot_ly(
  data = pca_df_3d,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~type,
  colors = c("#1f77b4", "#ff7f0e"),
  marker = list(size = 4)
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )

cat("PCA analysis completed.\n")
cat("Variance explained:\n")
print(summary(pca)$importance[2, ])
