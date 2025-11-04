#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-10-30
# Description: 

# ==== Main Script Start ====
#!/usr/bin/env Rscript
# Compute theta by matching edges via MRCA, then annotate tree
# Usage: Rscript annotate_theta_tree.R astral.tre iqtree.tre output.pdf

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
    library(ggplot2)
  library(dplyr)
  library(viridis)
    library(viridisLite)
})

args <- commandArgs(TRUE)
if(length(args) < 3){
  stop("Usage: Rscript annotate_theta_tree.R <ASTRAL_tree> <IQTREE_tree> <output.pdf>")
}

astral_file <- args[1]
iqtree_file <- args[2]
pdf_out <- args[3]

cat("Reading trees...\n")
tree_coal <- read.tree(astral_file)  # coalescent units (ASTRAL)
tree_mut  <- read.tree(iqtree_file)   # mutation units (IQTREE)

# ========= Match edges via MRCA =========

get_edge_label <- function(tree) {
  apply(tree$edge, 1, function(e) {
    tips <- sort(tree$tip.label[tree$edge[,2] <= length(tree$tip.label)])
    desc <- sort(tree$tip.label[getDescendants(tree, e[2], "tips")])
    paste(desc, collapse="|")
  })
}

getDescendants <- function(tree, node, type = c("tips","all")){
  type <- match.arg(type)
  kids <- tree$edge[tree$edge[,1]==node,2]
  if(length(kids)==0) return(NULL)
  desc <- kids
  for(k in kids) desc <- c(desc, getDescendants(tree, k, type))
  if(type=="tips") desc <- desc[desc <= length(tree$tip.label)]
  desc
}

coal_key <- get_edge_label(tree_coal)
mut_key  <- get_edge_label(tree_mut)

match_idx <- match(mut_key, coal_key)

if(any(is.na(match_idx))){
  stop("Cannot match edges between trees (likely tree topology mismatch)")
}

theta <- tree_mut$edge.length / tree_coal$edge.length[match_idx]

df <- data.frame(
  edge = 1:nrow(tree_mut$edge),
  parent = tree_mut$edge[,1],
  child = tree_mut$edge[,2],
  mutation_unit = tree_mut$edge.length,
  coalescent_unit = tree_coal$edge.length[match_idx],
  theta = theta
)

Ntip <- length(tree_mut$tip.label)
df_internal <- df[df$child > Ntip, ]
write.table(df_internal, "theta_per_node.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat("Theta file saved: theta_per_node.tsv\n")

# ========= Plot =========

p <- ggtree(tree_mut)
df_tree <- p$data %>% left_join(df_internal, by=c("node"="child"))

pdf(pdf_out, width=8, height=12)

p <- ggtree(tree_mut) %<+% df_tree +
    geom_tiplab(size=3, hjust=0) +
  geom_tree(aes(color = theta), size = 1) +
  scale_color_viridis(option="plasma", na.value="grey70") +
  geom_text2(aes(label = ifelse(!is.na(theta), round(theta, 3), "")),
             hjust=-0.2, size=2.5) +
  theme_tree2() +
  ggtitle("Theta mapped to tree")

print(p)
dev.off()

cat("Annotated tree saved to:", pdf_out, "\n")

