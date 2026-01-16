#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-08-27
# Description: Extract ASTRAL support values with node IDs
# usage: Rscript Get_support_value.R input_tree [output_csv]
# ==== Main Script Start ====

suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Get_support_value.R input_tree [output_csv]")
}

input_file <- args[1]
output_file <- if (length(args) >= 2) args[2] else "Support_values.csv"

if (!file.exists(input_file)) {
  stop("Error: Input file '", input_file, "' not found")
}

tree <- read.tree(input_file)

n_tips <- Ntip(tree)
node_ids <- (n_tips + 1):(n_tips + tree$Nnode)

df <- data.frame(
  Node = node_ids,
  Support = tree$node.label,
  stringsAsFactors = FALSE
)

df <- df[!is.na(df$Support) & df$Support != "", ]

write.csv(df, output_file, row.names = FALSE, quote = FALSE)

cat("Successfully exported", nrow(df), "support values to", output_file, "\n")