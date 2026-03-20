# !/usr/bin/Rscript
# Author: XiongTao
# Description: Compare RF distances between chloroplast tree and different factor-based nuclear gene trees

library(ape)
library(dplyr)
library(ggplot2)

# Load functions
source("utils_functions.R")

# Read chloroplast tree
chloroplast <- read.tree("tree/rosa_chloroplast_partition_rt.tre")

# Define factors to compare
factors <- c("Random_nuclear", "GC_content", "Mean_evolutionary_rate_without0", "PC1",
             "Proportion_parsimony_informative", "Proportion_variable_sites", "RCV",
             "Mean_internal_branch", "Mean_teminal_branch")

# Calculate RF distances for each factor
df_list <- list()

for (factor_name in factors) {
  tree_path <- paste0("./tree/", factor_name, "/")
  if (dir.exists(tree_path)) {
    tree_list <- list.files(tree_path, pattern = "\\.tre$", full.names = T)
    if (length(tree_list) > 0) {
      df_list[[factor_name]] <- get_rf_distanch(tree_list, chloroplast, factor_name)
      cat("Processed:", factor_name, "- Trees:", length(tree_list), "\n")
    }
  }
}

# Combine all results
df <- bind_rows(df_list)
df$type <- factor(df$type, levels = factors)

# Calculate baseline (median of random nuclear)
med <- median(df_list[["Random_nuclear"]]$rf)

# Create boxplot
p <- ggplot(df, aes(x = type, y = rf, fill = type)) +
  geom_boxplot() +
  geom_hline(yintercept = med, color = "red", linetype = "dashed") +
  scale_fill_manual(values = c("Random_nuclear" = "grey90", "GC_content" = "grey70",
                                "Mean_evolutionary_rate_without0" = "grey70", "PC1" = "grey70",
                                "Proportion_parsimony_informative" = "grey70", "Proportion_variable_sites" = "grey70",
                                "RCV" = "grey70", "Mean_internal_branch" = "grey70", "Mean_teminal_branch" = "grey70")) +
  labs(x = "", y = "RF Distance to Chloroplast Tree", title = "Comparison of RF Distances Across Different Factors") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
ggsave("Boxplot_all_diff.pdf", p, width = 20, height = 10)

cat("\nFinal comparison completed.\n")
cat("Baseline (Random nuclear median RF):", round(med, 4), "\n")
