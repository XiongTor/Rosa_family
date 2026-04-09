# !/usr/bin/Rscript
# Author: XiongTao
# Description: Identify nuclear genes most similar to chloroplast tree and visualize their distribution

library(ape)
library(phangorn)
library(dplyr)
library(ggplot2)
library(cowplot)

# Load functions and data
source("00-utils_functions.R")
merged_orthofinder <- readRDS("merged_orthofinder.rds")
merged_chloroplast <- readRDS("merged_chloroplast.rds")

# Read chloroplast tree as reference
chloroplast_tree <- read.tree("./tree/rosa_chloroplast_partition_rt.tre")

# Get all nuclear gene tree files
all_tree <- list.files("./tree/orthofinder_genes/", pattern = "\\.tre$", full.names = TRUE)

cat("Total nuclear gene trees:", length(all_tree), "\n")

# Calculate RF distance between each nuclear gene tree and chloroplast tree
rf_results <- data.frame(
  tree_file = character(),
  gene_name = character(),
  rf_distance = numeric(),
  stringsAsFactors = FALSE
)

cat("Calculating RF distances...\n")
for (i in seq_along(all_tree)) {
  if (i %% 100 == 0) cat("  Processed", i, "/", length(all_tree), "\n")

  tree_file <- all_tree[i]
  gene_name <- sub("_.*$", "", basename(tree_file))

  tryCatch({
    nuc_tree <- read.tree(tree_file)

    # Find common tips
    common_tips <- intersect(nuc_tree$tip.label, chloroplast_tree$tip.label)

    if (length(common_tips) >= 4) {
      # Prune trees to common tips
      nuc_pruned <- drop.tip(nuc_tree, setdiff(nuc_tree$tip.label, common_tips))
      chloro_pruned <- drop.tip(chloroplast_tree, setdiff(chloroplast_tree$tip.label, common_tips))

      # Calculate RF distance
      rf_dist <- RF.dist(nuc_pruned, chloro_pruned, normalize = TRUE)

      rf_results <- rbind(rf_results, data.frame(
        tree_file = tree_file,
        gene_name = gene_name,
        rf_distance = rf_dist,
        stringsAsFactors = FALSE
      ))
    }
  }, error = function(e) {
    warning(paste("Error processing", gene_name, ":", e$message))
  })
}

cat("RF distance calculation completed.\n")
cat("Valid comparisons:", nrow(rf_results), "\n")

# Select top 20% genes with smallest RF distance
threshold <- quantile(rf_results$rf_distance, 0.2, na.rm = TRUE)
top20_genes <- rf_results %>%
  filter(rf_distance <= threshold) %>%
  arrange(rf_distance)

cat("\nTop 20% threshold RF distance:", threshold, "\n")
cat("Number of genes in top 20%:", nrow(top20_genes), "\n")

# Merge with merged_orthofinder to get attribute values
data_like <- merged_orthofinder %>%
  filter(Alignment_name %in% top20_genes$gene_name)

cat("Genes found in merged_orthofinder:", nrow(data_like), "\n")

# Select bottom 20% genes with largest RF distance (most unlike chloroplast)
threshold_unlike <- quantile(rf_results$rf_distance, 0.8, na.rm = TRUE)
bottom20_genes <- rf_results %>%
  filter(rf_distance >= threshold_unlike) %>%
  arrange(desc(rf_distance))

cat("\nBottom 20% threshold RF distance:", threshold_unlike, "\n")
cat("Number of genes in bottom 20%:", nrow(bottom20_genes), "\n")

# Merge with merged_orthofinder to get attribute values
data_unlike <- merged_orthofinder %>%
  filter(Alignment_name %in% bottom20_genes$gene_name)

cat("Genes found in merged_orthofinder (unlike):", nrow(data_unlike), "\n")

# Define attributes to plot
attributes <- c(
  "GC_content",
  "Proportion_variable_sites",
  "Proportion_parsimony_informative",
  "Mean_evolutionary_rate",
  "Mean_evolutionary_rate_without0",
  "RCV",
  "Mean_internal_branch",
  "Mean_teminal_branch"
)

# Create ranked bar plots for each attribute
plot_list_bar <- list()

for (attr in attributes) {
  cat("Plotting ranked bar for", attr, "...\n")

  # Sort merged_orthofinder by attribute
  merged_sorted <- merged_orthofinder %>%
    arrange(.data[[attr]]) %>%
    mutate(
      Rank = row_number(),
      Group = ifelse(Alignment_name %in% data_like$Alignment_name,
                     "Chloroplast-like",
                     "Others")
    )

  # Create plot
  p <- ggplot(merged_sorted, aes(x = Rank, y = .data[[attr]], fill = Group)) +
    geom_col(width = 1) +
    scale_fill_manual(
      values = c("Chloroplast-like" = "#E41A1C", "Others" = "grey80"),
      breaks = c("Chloroplast-like", "Others")
    ) +
    labs(
      title = attr,
      x = "Genes (Ranked)",
      y = attr
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )

  plot_list_bar[[attr]] <- p
}

# Add legend to the last plot
plot_list_bar[[length(plot_list_bar)]] <- plot_list_bar[[length(plot_list_bar)]] +
  theme(legend.position = "bottom")

# Create density plots for each attribute
plot_list_density <- list()

for (attr in attributes) {
  cat("Plotting density for", attr, "...\n")

  # Get x-axis range
  all_values <- c(
    merged_orthofinder[[attr]],
    merged_chloroplast[[attr]],
    data_like[[attr]],
    data_unlike[[attr]]
  )
  x_min <- min(all_values, na.rm = TRUE)
  x_max <- max(all_values, na.rm = TRUE)

  # Create density plot
  p <- ggplot() +
    geom_density(data = merged_orthofinder, aes(x = .data[[attr]], fill = "All nuclear"),
                 adjust = 3, alpha = 0.4) +
    geom_density(data = merged_chloroplast, aes(x = .data[[attr]], fill = "Chloroplast"),
                 adjust = 3, alpha = 0.6) +
    geom_density(data = data_like, aes(x = .data[[attr]], fill = "Chloroplast-like"),
                 adjust = 3, alpha = 0.6) +
    geom_density(data = data_unlike, aes(x = .data[[attr]], fill = "Chloroplast-unlike"),
                 adjust = 3, alpha = 0.6) +
    scale_fill_manual(
      name = "",
      values = c(
        "All nuclear" = "yellow",
        "Chloroplast" = "#a64d6d",
        "Chloroplast-like" = "#2073a5",
        "Chloroplast-unlike" = "#4DAF4A"
      ),
      breaks = c("Chloroplast", "All nuclear", "Chloroplast-like", "Chloroplast-unlike")
    ) +
    labs(x = attr, y = "Density", title = attr) +
    xlim(x_min, x_max) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )

  plot_list_density[[attr]] <- p
}

# Add legend to the last density plot
plot_list_density[[length(plot_list_density)]] <- plot_list_density[[length(plot_list_density)]] +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# Combine all bar plots in 2 columns
final_plot_bar <- plot_grid(
  plotlist = plot_list_bar,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# Combine all density plots in 2 columns
final_plot_density <- plot_grid(
  plotlist = plot_list_density,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

# Save plots
pdf("chloroplast_like_genes_ranked_distribution.pdf", width = 12, height = 16)
print(final_plot_bar)
dev.off()

pdf("chloroplast_like_genes_density_distribution.pdf", width = 12, height = 16)
print(final_plot_density)
dev.off()

cat("\nPlots saved:\n")
cat("  - chloroplast_like_genes_ranked_distribution.pdf\n")
cat("  - chloroplast_like_genes_density_distribution.pdf\n")

# Save results
write.csv(top20_genes, "top20_chloroplast_like_genes.csv", row.names = FALSE, quote = FALSE)
write.csv(data_like, "chloroplast_like_genes_attributes.csv", row.names = FALSE, quote = FALSE)
write.csv(bottom20_genes, "bottom20_chloroplast_unlike_genes.csv", row.names = FALSE, quote = FALSE)
write.csv(data_unlike, "chloroplast_unlike_genes_attributes.csv", row.names = FALSE, quote = FALSE)

cat("\nResults saved:\n")
cat("  - top20_chloroplast_like_genes.csv\n")
cat("  - chloroplast_like_genes_attributes.csv\n")
cat("  - bottom20_chloroplast_unlike_genes.csv\n")
cat("  - chloroplast_unlike_genes_attributes.csv\n")
cat("\nAnalysis completed!\n")
