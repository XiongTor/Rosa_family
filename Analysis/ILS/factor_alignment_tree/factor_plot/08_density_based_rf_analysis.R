# !/usr/bin/Rscript
# Author: XiongTao
# Description: RF distance analysis using density-based equidistant intervals

library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(cowplot)

# Load functions and data
source("00-utils_functions.R")
merged_orthofinder <- readRDS("merged_orthofinder.rds")

# Set tree path
all_tree <- list.files("./tree/orthofinder_genes/", pattern = "\\.tre$", full.names = T)

# Configuration
factor <- "Mean_internal_branch"  # Change this to analyze different factors
interval_mode <- "paired"  # "paired" or "consecutive"

# Define quantile points and intervals
# "paired": 隔一个取一段 (1-2, 3-4, 5-6, ...)
# "consecutive": 取所有连续段 (1-2, 2-3, 3-4, 4-5, ...)
quantile_points <- c(0.40, 0.45, 0.5, 0.55, 0.605, 0.65, 0.71, 0.76, 0.82, 0.86)
n_intervals <- 5
interval_colors <- c("grey", "#4DAF4A", "#377EB8", "orange", "#E41A1C")
interval_labels <- c("Low", "Low_m", "Mid", "Mid_h", "High")

# Prepare data
merged <- merged_orthofinder

# Calculate density
d <- density(merged[[factor]], adjust = 3)
df_density <- data.frame(x = d$x, y = d$y)

# Calculate actual quantile values based on data range
x <- merged[[factor]]
qs <- round(
  min(x, na.rm = TRUE) +
    quantile_points * (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
  2
)

# Extract gene lists for each interval based on density distribution
interval_data <- list()
genes_list <- list()

for (i in 1:n_intervals) {
  # Determine start and end indices based on mode
  if (interval_mode == "paired") {
    idx_start <- 2 * i - 1
    idx_end <- 2 * i
  } else {  # consecutive
    idx_start <- i
    idx_end <- i + 1
  }

  # Extract interval data for density plot
  interval_data[[i]] <- df_density[df_density$x >= qs[idx_start] & df_density$x <= qs[idx_end], ]

  # Extract gene names falling within this value range
  genes_list[[i]] <- merged[["Alignment_name"]][
    merged[[factor]] >= qs[idx_start] &
      merged[[factor]] <= qs[idx_end]
  ]

  # Set names
  if (length(interval_labels) >= i) {
    names(genes_list)[i] <- interval_labels[i]
  } else {
    names(genes_list)[i] <- paste0("Interval_", i)
  }
}

# Print summary
cat("\nDensity-based interval extraction summary:\n")
for (i in 1:n_intervals) {
  cat(sprintf("  %s: %d genes (value range: %.4f - %.4f)\n",
              names(genes_list)[i],
              length(genes_list[[i]]),
              qs[2*i-1],
              qs[2*i]))
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

# Build factor value dataframe for violin plot
df_factor <- do.call(rbind, lapply(1:n_intervals, function(i) {
  data.frame(
    Factor_value = merged[[factor]][
      merged[["Alignment_name"]] %in% genes_list[[i]]
    ],
    Group = names(genes_list)[i]
  )
}))

df_factor$Group <- factor(df_factor$Group, levels = names(genes_list))

# Create plots

# Plot 1: Density distribution with highlighted intervals
p1 <- ggplot(df_density, aes(x, y)) +
  geom_area(fill = "#7393B3", alpha = 0.75) +
  geom_line(linewidth = 0.6) +
  theme_bw() +
  labs(x = factor, y = "Density", title = paste("Density-based partitioning:", factor))

# Add interval ribbons and labels
for (i in 1:n_intervals) {
  idx_start <- if (interval_mode == "paired") 2 * i - 1 else i
  idx_end <- if (interval_mode == "paired") 2 * i else i + 1

  p1 <- p1 +
    geom_ribbon(
      data = interval_data[[i]],
      aes(ymin = 0, ymax = y),
      fill = interval_colors[i],
      alpha = 0.6
    ) +
    annotate("text", x = qs[idx_start], y = 0,
             label = qs[idx_start],
             vjust = 1.5, size = 3) +
    annotate("text", x = qs[idx_end], y = 0,
             label = qs[idx_end],
             vjust = 1.5, size = 3)
}

# Plot 2: Factor value distribution by group
p2 <- ggplot(df_factor, aes(x = Group, y = Factor_value, fill = Group)) +
  geom_violin(width = 1, alpha = 0.4, trim = TRUE) +
  geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
  scale_fill_manual(values = interval_colors) +
  theme_classic(base_size = 14) +
  ylab(factor) + xlab("") +
  theme(legend.position = "none")

# Plot 3: RF distance distribution by group
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

# Combine and save plots
pdf(paste0(factor, "_density_based_equidistants.pdf"), width = 12, height = 10)
plot_grid(
  p1,
  plot_grid(p2, p3, ncol = 2, align = "h", axis = "tb", rel_widths = c(1, 1)),
  ncol = 1, rel_heights = c(1, 1)
)
dev.off()

cat("\nDensity-based RF distance analysis completed for factor:", factor, "\n")
cat("Output saved to:", paste0(factor, "_density_based_equidistants.pdf\n"))
