# !/usr/bin/Rscript
# Author: XiongTao
# Description: Correlation analysis and visualization

library(dplyr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(ggplot2)

# Load data
merged_orthofinder <- readRDS("merged_orthofinder.rds")

# Define all columns and differential factors
cols <- c(
  "No_of_taxa", "Alignment_length", "GC_content",
  "Proportion_variable_sites", "Proportion_parsimony_informative",
  "Mean_evolutionary_rate", "Mean_evolutionary_rate_without0",
  "RCV", "Mean_bipartition", "Mean_internal_branch",
  "Mean_teminal_branch", "treeness", "saturation",
  "absolute_saturation", "treeness.RCV"
)

diff <- c("GC_content", "Proportion_variable_sites",
          "Proportion_parsimony_informative", "Mean_evolutionary_rate", "RCV")

# Calculate correlation
df <- merged_orthofinder[, cols]
res <- rcorr(as.matrix(df))

cor_mat <- res$r
p_mat <- res$P

# Plot correlation matrix
pdf("corrplot_diff_all_factor.pdf", width = 10, height = 10)
corrplot(
  cor_mat,
  method = "square",
  tl.col = "black",
  p.mat = p_mat,
  insig = "blank",
  addCoef.col = "black",
  number.cex = 1
)
dev.off()

# Create column plot for differential factors
df <- merged_orthofinder %>%
  arrange(desc(Proportion_variable_sites)) %>%
  mutate(order = row_number()) %>%
  select(order, all_of(diff))

col <- ggplot(df, aes(x = order)) +
  geom_col(aes(y = Proportion_variable_sites, fill = "PVS"), alpha = 0.5) +
  geom_col(aes(y = Proportion_parsimony_informative, fill = "PPI"), alpha = 0.6) +
  geom_col(aes(y = Mean_evolutionary_rate, fill = "Rate"), alpha = 1) +
  geom_col(aes(y = GC_content, fill = "GC"), alpha = 0.6) +
  scale_fill_manual(
    breaks = c("PVS", "PPI", "Rate", "GC"),
    values = c("PVS" = "#0072B2", "PPI" = "red", "Rate" = "black", "GC" = "#DAA520"),
    labels = c("PVS" = "Variable sites", "PPI" = "Parsimony informative",
               "Rate" = "Evolutionary rate", "GC" = "GC content")
  ) +
  labs(x = "Genes", y = "Value") +
  theme_bw() +
  theme(
    legend.position = c(0.88, 0.85),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"),
    legend.title = element_blank()
  )

ggsave("col.pdf", col, width = 7, height = 5)

cat("Correlation analysis completed.\n")
