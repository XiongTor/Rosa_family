# !/usr/bin/Rscript
# Author: XiongTao
# date: 2025.12.11
# Description: Utility functions for Rosaceae cytonuclear analysis

# Function: Calculate pairwise RF distance for one column
pairwise_rf_onecol <- function(
    genes,
    all_tree,
    min_tips = 4,
    pattern_suffix = "_",
    normalize = TRUE,
    rooted = FALSE
) {
  require(ape)
  require(phangorn)

  matched_files <- all_tree[
    sapply(all_tree, function(f)
      any(startsWith(basename(f), paste0(genes, pattern_suffix)))
    )
  ]

  if (length(matched_files) < 2) {
    warning("Matched tree files < 2, RF not computed.")
    return(matrix(NA, ncol = 1, dimnames = list(NULL, "RF")))
  }

  gene_trees <- lapply(matched_files, read.tree)
  names(gene_trees) <- sub("_.*$", "", basename(matched_files))

  n <- length(gene_trees)
  rf_vec <- c()

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      common_tips <- intersect(
        gene_trees[[i]]$tip.label,
        gene_trees[[j]]$tip.label
      )

      if (length(common_tips) < min_tips) next

      ti <- drop.tip(
        gene_trees[[i]],
        setdiff(gene_trees[[i]]$tip.label, common_tips)
      )
      tj <- drop.tip(
        gene_trees[[j]],
        setdiff(gene_trees[[j]]$tip.label, common_tips)
      )

      rf_vec <- c(
        rf_vec,
        RF.dist(
          ti, tj,
          normalize = normalize,
          rooted = rooted
        )
      )
    }
  }

  rf_mat <- matrix(rf_vec, ncol = 1)
  colnames(rf_mat) <- "RF"

  return(rf_mat)
}

# Function: Match chloroplast property once
match_property_once <- function(chloro, nuc, prop, pct=0.25){

  chloro <- chloro[order(chloro[[prop]]), ]
  nuc    <- nuc[order(nuc[[prop]]), ]

  selected <- c()

  result <- data.frame(
    chloro_gene   = character(),
    chloro_value  = numeric(),
    matched_gene  = character(),
    matched_value = numeric(),
    diff          = numeric(),
    stringsAsFactors = FALSE
  )

  for(i in seq_len(nrow(chloro))){

    target_val <- chloro[[prop]][i]
    lower <- target_val  - abs(target_val*pct)
    upper <- target_val + abs(target_val*pct)

    cand_idx <- which(
      nuc[[prop]] >= lower &
        nuc[[prop]] <= upper &
        !(nuc$Alignment_name %in% selected)
    )

    if(length(cand_idx)==0){
      matched_gene  <- NA
      matched_value <- NA
      diff <- NA

    } else {

      diffs <- abs(nuc[[prop]][cand_idx] - target_val)
      picked_idx <- cand_idx[which.min(diffs)]

      matched_gene  <- nuc$Alignment_name[picked_idx]
      matched_value <- nuc[[prop]][picked_idx]
      diff <- abs(target_val - matched_value)
      selected <- c(selected, matched_gene)
    }

    result <- rbind(
      result,
      data.frame(
        chloro_gene   = chloro$Alignment_name[i],
        chloro_value  = target_val,
        matched_gene  = matched_gene,
        matched_value = matched_value,
        diff          = diff
      )
    )
  }

  return(result)
}

# Function: Match property once with random sampling capability
match_property_once_random <- function(chloro, nuc, prop, pct=0.25){

  chloro <- chloro[order(chloro[[prop]]), ]
  nuc    <- nuc[order(nuc[[prop]]), ]

  result <- data.frame(
    chloro_gene   = character(),
    chloro_value  = numeric(),
    matched_gene  = character(),
    matched_value = numeric(),
    diff          = numeric(),
    all_nuc_gene  = character(),
    n_matched     = integer(),
    stringsAsFactors = FALSE
  )

  selected <- c()
  first_match <- list()

  for(i in seq_len(nrow(chloro))){
    target_val <- chloro[[prop]][i]
    lower <- target_val  - abs(target_val*pct)
    upper <- target_val + abs(target_val*pct)

    cand_idx <- which(
      nuc[[prop]] >= lower &
        nuc[[prop]] <= upper &
        !(seq_len(nrow(nuc)) %in% selected)
    )

    if(length(cand_idx) == 0){
      first_match[[i]] <- list(idx = NULL, gene = NA, value = NA, diff = NA)
    } else {
      diffs <- abs(target_val - nuc[[prop]][cand_idx])
      best_idx <- cand_idx[which.min(diffs)]

      first_match[[i]] <- list(
        idx   = best_idx,
        gene  = nuc$Alignment_name[best_idx],
        value = nuc[[prop]][best_idx],
        diff  = abs(target_val - nuc[[prop]][best_idx])
      )

      selected <- c(selected, best_idx)
    }
  }

  second_match <- list()

  for(i in seq_len(nrow(chloro))){
    target_val <- chloro[[prop]][i]
    lower <- target_val  - abs(target_val*pct)
    upper <- target_val + abs(target_val*pct)

    cand_idx <- which(
      nuc[[prop]] >= lower &
        nuc[[prop]] <= upper &
        !(seq_len(nrow(nuc)) %in% selected)
    )

    if(length(cand_idx) == 0){
      second_match[[i]] <- list(idx = NULL, gene = NA)
    } else {
      diffs <- abs(target_val - nuc[[prop]][cand_idx])
      best_idx <- cand_idx[which.min(diffs)]

      second_match[[i]] <- list(
        idx  = best_idx,
        gene = nuc$Alignment_name[best_idx]
      )

      selected <- c(selected, best_idx)
    }
  }

  for(i in seq_len(nrow(chloro))){
    first  <- first_match[[i]]
    second <- second_match[[i]]

    all_genes <- c()
    if(!is.na(first$gene)) all_genes <- c(all_genes, first$gene)
    if(!is.na(second$gene)) all_genes <- c(all_genes, second$gene)

    result <- rbind(
      result,
      data.frame(
        chloro_gene   = chloro$Alignment_name[i],
        chloro_value  = chloro[[prop]][i],
        matched_gene  = first$gene,
        matched_value = first$value,
        diff          = first$diff,
        all_nuc_gene  = if(length(all_genes) > 0) paste(all_genes, collapse = "; ") else NA,
        n_matched     = length(all_genes)
      )
    )
  }

  return(result)
}

# Function: Random sample from candidates
random_sample_from_candidates <- function(candidates_df, n_iterations=20){

  all_results <- list()

  for(iter in 1:n_iterations){

    iter_result <- data.frame(
      chloro_gene   = character(),
      matched_gene  = character(),
      iteration     = integer(),
      stringsAsFactors = FALSE
    )

    for(i in seq_len(nrow(candidates_df))){

      chloro_gene <- candidates_df$chloro_gene[i]
      all_nuc <- candidates_df$all_nuc_gene[i]

      if(is.na(all_nuc)){
        matched_gene <- NA
      } else {
        nuc_genes <- strsplit(all_nuc, "; ")[[1]]
        matched_gene <- sample(nuc_genes, 1)
      }

      iter_result <- rbind(
        iter_result,
        data.frame(
          chloro_gene  = chloro_gene,
          matched_gene = matched_gene,
          iteration    = iter,
          stringsAsFactors = FALSE
        )
      )
    }

    all_results[[iter]] <- iter_result
  }

  final_result <- do.call(rbind, all_results)

  return(final_result)
}

# Function: Plot random iterations
plot_random_iterations <- function(random_results, nuc_data, chloro_data, all_nuc_data, prop = "Mean_evolutionary_rate"){

  p <- ggplot() +
    geom_density(data = all_nuc_data, aes(x = .data[[prop]]),
                 adjust = 3, fill = "yellow", alpha = 0.6) +
    geom_density(data = chloro_data, aes(x = .data[[prop]]),
                 adjust = 3, fill = "#a64d6d", alpha = 0.6)

  for(iter in unique(random_results$iteration)){

    iter_genes <- random_results %>%
      filter(iteration == iter, !is.na(matched_gene)) %>%
      pull(matched_gene)

    iter_data <- nuc_data %>%
      filter(Alignment_name %in% iter_genes)

    if(nrow(iter_data) > 0){
      p <- p + geom_density(data = iter_data, aes(x = .data[[prop]]),
                            adjust = 3, fill = "#2073a5", alpha = 0.1)
    }
  }

  p <- p +
    labs(x = prop, y = "Density",
         title = paste0("Distribution of ", prop, " across 20 random iterations")) +
    theme_bw()

  return(p)
}

# Function: Sample every n genes
sample_every_n_genes <- function(chloro_data, prop = "Mean_evolutionary_rate", interval = 3){

  library(dplyr)

  chloro_sorted <- chloro_data %>%
    arrange(.data[[prop]])

  n_total <- nrow(chloro_sorted)
  selected_indices <- seq(1, n_total, by = interval)

  result <- chloro_sorted[selected_indices, ]

  cat("Sampling summary:\n")
  cat(sprintf("  Total genes in dataset: %d\n", n_total))
  cat(sprintf("  Sampling interval: every %d genes\n", interval))
  cat(sprintf("  Sampled genes: %d\n", nrow(result)))
  cat(sprintf("  GC_content range: %.4f - %.4f\n",
              min(result[[prop]]), max(result[[prop]])))

  return(result)
}

# Function: Get RF distance between trees and chloroplast
get_rf_distanch <- function(tree_list, chloroplast, name){
  rf_list <- numeric(length(tree_list))
  for (i in seq_along(tree_list)){
    tt <- read.tree(tree_list[i])
    rf_list[i] <- RF.dist(tt, chloroplast, check.labels = TRUE, normalize = TRUE)
  }
  df_chloroplast <- data.frame(rf = rf_list)
  df_chloroplast$type <- name
  return(df_chloroplast)
}

# Function: Auto-tune pct value for optimal matching
auto_tune_pct <- function(chloro, nuc, prop,
                          target_n = NULL,
                          pct_start = 0.1,
                          pct_max = 3.0,
                          step = 0.05,
                          tolerance = 3) {

  if(is.null(target_n)) {
    target_n <- nrow(chloro)
  }

  cat(sprintf("Target: match %d genes (±%d)\n", target_n, tolerance))
  cat(sprintf("Property: %s\n", prop))
  cat("Searching for optimal pct value...\n\n")

  pct_values <- seq(pct_start, pct_max, by = step)
  results <- data.frame(
    pct = numeric(),
    n_matched = integer(),
    n_with_2_candidates = integer()
  )

  for(p in pct_values) {
    result <- match_property_once_random(chloro, nuc, prop, pct = p)

    n_matched <- sum(!is.na(result$matched_gene))
    n_with_2 <- sum(result$n_matched == 2, na.rm = TRUE)

    results <- rbind(results, data.frame(
      pct = p,
      n_matched = n_matched,
      n_with_2_candidates = n_with_2
    ))

    cat(sprintf("pct=%.2f: matched %d genes (%d with 2 candidates)\n",
                p, n_matched, n_with_2))

    if(n_matched >= target_n && n_matched <= target_n + tolerance) {
      cat(sprintf("\n✓ Optimal pct found: %.2f\n", p))
      cat(sprintf("  Matched genes: %d\n", n_matched))
      cat(sprintf("  Genes with 2 candidates: %d\n", n_with_2))
      return(list(optimal_pct = p, n_matched = n_matched, results = results))
    }
  }

  results$diff <- abs(results$n_matched - target_n)
  best_idx <- which.min(results$diff)
  best_pct <- results$pct[best_idx]

  cat(sprintf("\n⚠ No exact match found, returning closest: pct = %.2f\n", best_pct))
  cat(sprintf("  Matched genes: %d (target: %d)\n",
              results$n_matched[best_idx], target_n))

  return(list(
    optimal_pct = best_pct,
    n_matched = results$n_matched[best_idx],
    results = results
  ))
}
