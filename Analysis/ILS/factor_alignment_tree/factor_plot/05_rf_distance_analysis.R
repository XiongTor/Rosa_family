# !/usr/bin/Rscript
# Author: XiongTao
# Description: RF distance analysis for gene partitions

library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(cowplot)

# Load functions and data
source("utils_functions.R")
merged_orthofinder <- readRDS("merged_orthofinder.rds")

# Set tree path
all_tree <- list.files("./tree/orthofinder_genes/", pattern = "\\.tre$", full.names = T)

# Configuration
factor <- "Mean_internal_branch"  # Change this to analyze different factors
interval_mode <- "paired"  # "paired" or "consecutive"

# Define quantile points and intervals
quantile_points <- c(0.1, 0.15, 0.25, 0.3, 0.4, 0.45, 0.55, 0.6, 0.7, 0.75)
n_intervals <- 5
interval_colors <- c("grey", "#4DAF4A", "#377EB8", "orange", "#E41A1C")
interval_labels <- c("Low", "Low_m", "Mid", "Mid_h", "High")

# Sort data by factor
merged <- merged_orthofinder[order(merged_orthofinder[[factor]]), ]
n_total <- nrow(merged)

# Extract gene lists for each interval
genes_list <- list()
df_pvs <- data.frame()

for (group_name in interval_labels) {
  idx <- which(interval_labels == group_name)
  start <- quantile_points[2 * idx - 1]
  end <- quantile_points[2 * idx]

  start_index <- floor(n_total * start) + 1
  end_index <- floor(n_total * end)

  tmp_df <- merged[start_index:end_index, c("Alignment_name", factor)]
  tmp_df$Group <- group_name
  df_pvs <- rbind(df_pvs, tmp_df)
  genes_list[[group_name]] <- merged$Alignment_name[start_index:end_index]
}

# Calculate RF distances
rf_list <- lapply(names(genes_list), function(group_name) {
  message(paste("Calculating RF distance for:", group_name))
  pairwise_rf_onecol(
    genes = genes_list[[group_name]],
    pattern_suffix = "_",
    all_tree = all_tree,
    min_tips = 5
  )
})
names(rf_list) <- names(genes_list)

# Build RF dataframe
df_rf <- do.call(rbind, lapply(names(genes_list), function(group_name) {
  data.frame(RF = as.numeric(rf_list[[group_name]]), Group = group_name)
}))

df_rf$Group <- factor(df_rf$Group, levels = names(genes_list))
df_pvs$Group <- factor(df_pvs$Group, levels = names(genes_list))

# Create plots
merged$Group <- "Others"
for (group_name in names(genes_list)) {
  idx <- which(interval_labels == group_name)
  start_pct <- quantile_points[2 * idx - 1]
  end_pct <- quantile_points[2 * idx]
  start_idx <- floor(nrow(merged) * start_pct) + 1
  end_idx <- floor(nrow(merged) * end_pct)
  merged$Group[start_idx:end_idx] <- group_name
}

merged$Group <- factor(merged$Group, levels = c(names(genes_list), "Others"))

p1 <- ggplot(merged, aes(x = reorder(Alignment_name, !!sym(factor)), y = .data[[factor]], fill = Group)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(interval_colors, "grey"), breaks = names(genes_list)) +
  labs(title = paste("Data Partitioning based on", factor), x = "Genes (Ranked)", y = factor) +
  theme_minimal() +
  theme(legend.position = "top")

p2 <- ggplot(df_pvs, aes(x = Group, y = .data[[factor]], fill = Group)) +
  geom_violin(width = 1, alpha = 0.4, trim = TRUE) +
  geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = interval_colors) +
  theme_classic(base_size = 14) +
  ylab(factor) + xlab("") +
  theme(legend.position = "none")

gene_n <- sapply(genes_list, length)
group_labels <- paste0(names(gene_n), "\n(", gene_n, ")")
names(group_labels) <- names(gene_n)

p3 <- ggplot(df_rf, aes(x = Group, y = RF, fill = Group)) +
  geom_violin(width = 1, alpha = 0.4, trim = TRUE) +
  geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = interval_colors) +
  scale_x_discrete(labels = group_labels) +
  theme_classic(base_size = 14) +
  ylab("RF") + xlab("") +
  theme(legend.position = "none")

pdf(paste0(factor, "_equidistants.pdf"), width = 12, height = 10)
plot_grid(p1, plot_grid(p2, p3, ncol = 2, align = "h", axis = "tb", rel_widths = c(1, 1)),
          ncol = 1, rel_widths = c(1, 1))
dev.off()

cat("RF distance analysis completed for factor:", factor, "\n")
