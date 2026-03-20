# !/usr/bin/Rscript
# Author: XiongTao
# Description: Match chloroplast-like nuclear genes

library(dplyr)
library(tidyverse)
library(ggplot2)

# Load functions and data
source("./script/00-utils_functions.R")
merged_orthofinder <- readRDS("merged_orthofinder.rds")
merged_chloroplast <- readRDS("merged_chloroplast.rds")

# Configuration
factor <- "PC1"  # Change this to analyze different factors
n_iterations <- 20
sample_size <- 75

# Ensure factor is numeric
merged_chloroplast[[factor]] <- as.numeric(as.character(merged_chloroplast[[factor]]))
merged_orthofinder[[factor]] <- as.numeric(as.character(merged_orthofinder[[factor]]))

# Auto-tune pct value
optimal_result <- auto_tune_pct(
  chloro = merged_chloroplast,
  nuc = merged_orthofinder,
  prop = factor,
  target_n = sample_size,
  pct_start = 0.1,
  pct_max = 10.0,
  step = 0.05,
  tolerance = 3
)
pct_value <- optimal_result$optimal_pct
cat(sprintf("\nUsing optimal pct value: %.2f\n\n", pct_value))

# Match chloroplast-like genes
result_gene <- match_property_once_random(merged_chloroplast, merged_orthofinder, factor, pct = pct_value)
random_results <- random_sample_from_candidates(result_gene, n_iterations = n_iterations)

# Plot random iterations
p <- plot_random_iterations(
  random_results = random_results,
  nuc_data = merged_orthofinder,
  chloro_data = merged_chloroplast,
  all_nuc_data = merged_orthofinder,
  prop = factor
)

# Add random nuclear gene sampling
set.seed(123)
sampling_results <- list()
for (i in 1:n_iterations) {
  sampled <- merged_orthofinder[sample(nrow(merged_orthofinder), sample_size), ]
  sampled_with_iter <- sampled %>%
    select(Alignment_name, all_of(factor)) %>%
    mutate(Iteration = i)
  sampling_results[[i]] <- sampled_with_iter
  p <- p + geom_density(data = sampled, aes(x = .data[[factor]]),
                        adjust = 3, fill = "yellow", alpha = 0.1)
}

# Finalize plot
p <- p +
  labs(x = factor, y = "Density",
       title = paste0("Distribution of ", factor, " ", n_iterations, " random iterations (n=", sample_size, ")")) +
  theme_bw() +
  geom_point(aes(x = -Inf, y = -Inf, color = "Random nuclear genes"), size = 0) +
  geom_point(aes(x = -Inf, y = -Inf, color = "Chloroplast like gene"), size = 0) +
  geom_point(aes(x = -Inf, y = -Inf, color = "Chloroplast genes"), size = 0) +
  scale_color_manual(
    name = "Dataset",
    values = c("Random nuclear genes" = "yellow", "Chloroplast like gene" = "#2073a5", "Chloroplast genes" = "#a64d6d"),
    labels = c("Random nuclear genes" = "Random nuclear genes", "Chloroplast like gene" = "Chloroplast like gene", "Chloroplast genes" = "Chloroplast genes")
  ) +
  guides(color = guide_legend(override.aes = list(size = 5, shape = 15, alpha = 0.6))) +
  theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = "white", color = "black"), legend.title = element_text(face = "bold"))

print(p)
ggsave(paste0("./distribution_factor_plot_20/", factor, "_random_iterations.pdf"), p, width = 10, height = 6)

# Export results
all_sampling_results <- bind_rows(sampling_results)
all_sampling_results <- all_sampling_results %>% select(Alignment_name, Iteration, all_of(factor))

all_sampling_results_width <- all_sampling_results %>%
  group_by(Iteration) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = row_id, names_from = Iteration, values_from = Alignment_name, names_prefix = "Iteration_") %>%
  select(-row_id)

random_results_width <- random_results %>%
  group_by(iteration) %>%
  mutate(row_id = row_number()) %>%
  ungroup() %>%
  pivot_wider(id_cols = row_id, names_from = iteration, values_from = matched_gene, names_prefix = "Iteration_") %>%
  select(-row_id)

write.csv(random_results, paste0("./distribution_factor_plot_20/", factor, "_Chloroplast_like_nuclear_gene.csv"), row.names = F, quote = F)
write.csv(random_results_width, paste0("./distribution_factor_plot_20/", factor, "_Chloroplast_like_nuclear_gene_genelist.csv"), row.names = F, quote = F)

cat("Chloroplast-like gene matching completed for factor:", factor, "\n")

