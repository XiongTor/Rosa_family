# !/usr/bin/Rscript
# Author: XiongTao
# Date: 2025.12.11
# Description: Master script to run all analyses sequentially

cat("=== Rosaceae Cytonuclear Factor Analysis Pipeline ===\n\n")

# Step 1: Data Loading
cat("Step 1: Loading and preprocessing data...\n")
source("01_data_loading.R")
cat("✓ Data loading completed\n\n")

# Step 2: Density Plots
cat("Step 2: Generating density plots...\n")
source("02_density_plots.R")
cat("✓ Density plots completed\n\n")

# Step 3: Correlation Analysis
cat("Step 3: Performing correlation analysis...\n")
source("03_correlation_analysis.R")
cat("✓ Correlation analysis completed\n\n")

# Step 4: PCA Analysis
cat("Step 4: Running PCA analysis...\n")
source("04_pca_analysis.R")
cat("✓ PCA analysis completed\n\n")

# Step 5: RF Distance Analysis - Quantile-based (optional - requires manual factor configuration)
cat("Step 5: RF distance analysis - Quantile-based (skipped - run manually with desired factor)\n")
cat("   Edit 05_rf_distance_analysis.R to set factor and parameters\n")
cat("   This uses equal gene count intervals (4.2 strategy)\n\n")

# Step 6: Chloroplast Matching (optional - requires manual factor configuration)
cat("Step 6: Chloroplast-like gene matching (skipped - run manually with desired factor)\n")
cat("   Edit 06_chloroplast_matching.R to set factor and pct value\n\n")

# Step 7: Final Comparison (optional - requires tree files)
cat("Step 7: Final comparison (skipped - run manually after tree construction)\n")
cat("   Run 07_final_comparison.R after building trees\n\n")

# Step 8: RF Distance Analysis - Density-based (optional - requires manual factor configuration)
cat("Step 8: RF distance analysis - Density-based (skipped - run manually with desired factor)\n")
cat("   Edit 08_density_based_rf_analysis.R to set factor and parameters\n")
cat("   This uses density distribution intervals (4.1 strategy)\n\n")

cat("=== Pipeline completed ===\n")
cat("Note: Steps 5-7 require manual configuration and tree files.\n")
