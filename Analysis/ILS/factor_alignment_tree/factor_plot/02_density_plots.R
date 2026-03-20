# !/usr/bin/Rscript
# Author: XiongTao
# Description: Generate density plots for all factors

library(ggplot2)
library(cowplot)

# Load data
merged_orthofinder <- readRDS("merged_orthofinder.rds")
merged_chloroplast <- readRDS("merged_chloroplast.rds")

# Define columns to plot
cols <- c(
  "No_of_taxa", "Alignment_length", "GC_content",
  "Proportion_variable_sites", "Proportion_parsimony_informative",
  "Mean_evolutionary_rate", "Mean_evolutionary_rate_without0",
  "RCV", "Mean_bipartition", "Mean_internal_branch",
  "Mean_teminal_branch", "treeness", "saturation",
  "absolute_saturation", "treeness.RCV"
)

# Create plots
plot_list <- list()
n_col <- 3

for (i in seq_along(cols)) {
  col_name <- cols[i]

  all_values <- c(merged_orthofinder[[col_name]], merged_chloroplast[[col_name]])
  x_min <- min(all_values, na.rm = TRUE)
  x_max <- max(all_values, na.rm = TRUE)

  p <- ggplot() +
    geom_density(data = merged_orthofinder, aes_string(x = col_name),
                 adjust = 3, fill = "#2073a5", alpha = 0.4) +
    geom_density(data = merged_chloroplast, aes_string(x = col_name),
                 adjust = 3, fill = "#a64d6d", alpha = 0.4) +
    labs(x = col_name, y = "Density") +
    xlim(x_min, x_max) +
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(size = 18),
      axis.text.y  = element_text(size = 18)
    )

  if (i == 1) {
    dens_max <- max(density(merged_orthofinder[[col_name]], na.rm = TRUE)$y,
                    density(merged_chloroplast[[col_name]], na.rm = TRUE)$y)
    p <- p +
      annotate("rect", xmin = x_min, xmax = x_min + (x_max - x_min) * 0.08,
               ymin = dens_max * 0.85, ymax = dens_max * 0.95, fill = "#2073a5") +
      annotate("rect", xmin = x_min, xmax = x_min + (x_max - x_min) * 0.08,
               ymin = dens_max * 0.65, ymax = dens_max * 0.75, fill = "#a64d6d") +
      annotate("text", x = x_min + (x_max - x_min) * 0.1,
               y = dens_max * 0.9, label = "orthofinder", hjust = 0, size = 5) +
      annotate("text", x = x_min + (x_max - x_min) * 0.1,
               y = dens_max * 0.7, label = "Chloroplast", hjust = 0, size = 5)
  } else {
    p <- p + theme(axis.title.y = element_blank())
  }

  plot_list[[col_name]] <- p
}

# Combine and save
final <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
ggsave("all_density_plots.pdf", final, width = 15, height = 15)

cat("Density plots saved to all_density_plots.pdf\n")
