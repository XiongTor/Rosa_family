#!/bin/bash
# Author: Tao Xiong
# Date: 2025-06-25
# Description: This script is used to analyze the relative contributions of ILS, GTEE, and IH to phylogenetic discordance.
# ==== 主体代码开始 ====

################################ 1. building phyloglenetic tree ###################################
# we need to build a species tree and many bootstrap replicates of gene tree. And all of them should be rooted.
# Get 100 bootstrap gene trees, it means we need 353*100=35300 gene trees.
for seq in `ls ../trimal/*.fasta`; do
  iqtree -s $seq -m MFP --bnni -T 40 -b 100 -pre ./$(basename $seq .fasta)
done

# Delete the branch length messages in speciestree
nw_topology -I rosa_Astral_species_rt.tre > rosa_Astral_species_rt_nobranch.tre

## And you should cat all genetrees into one file, and so does the bootstrap trees
cat *.genetrees > rosa_gene_trees.tre
cat *.bootstrapped.tre > rosa_gene_trees_bootstrap.tre


################################ 2. calculate the relative contributions of ILS, GTEE, and IH ##########

# Simulate gene trees under coalescent model for each BP species tree
# variations in the resulting simulated gene trees will reflect estimation error, ILS, and missing data in the empirical data set

mkdir geneTr_sim
Rscript --vanilla /data/xiongtao/tree/ILS/gene_flow/gene_flow_analy/script/MSC_geneTr_simulator.R $speciesTr $speciesTr_BP $geneTr
rm geneTr_sim/*.tem.genetrees
#excepted result: sets of simulated gene trees in the folder geneTr_sim


# Count triplet frequency in empirical gene trees and simulated gene trees
# This can take hours if there are >30 species, so it is strongly advised to distribute the work to the cluster. For example, submit a job to count triplet frequency for one set of gene trees.
# For empirical gene trees
python /data/xiongtao/tree/ILS/gene_flow/gene_flow_analy/script/triple_frequency_counter.py $geneTr $speciesTr
#excepted result: *.trp.tsv
#format: column1--species names of the triplet, sorted alphabetically; column2--triplet frequencies of (sp1,sp2);column3--triplet frequencies of (sp1,sp3);column4--triplet frequencies of (sp2,sp3).

# For simulated gene trees
# If the server hace enough cpu source, you can use parallel to speed up the process.

# !!!! Note, do not use the "./BP.sim.genetrees" as the input, it will make the output file become ".trp.tsv", just use "BP.sim.genetrees" is ok.
cat geneTr_sim/*.sim.genetrees>BP_sim.genetrees.txt
cat BP_sim.genetrees.txt |parallel -j 30 'python ./script/triple_frequency_counter.py {} rosa_Astral_species_rt.tre'



# Find significantly unbalanced triplets and map to species tree
# !!! The GitHub orgin script is wrong, the input file should be the genetrees,not the species tree. You can check the `find_unbalanced_triplets.py`
python /data/xiongtao/tree/ILS/gene_flow/gene_flow_analy/script/find_unbalanced_triplets.py $geneTr

#result in unbalanced.trp.tsv
python /data/xiongtao/tree/ILS/gene_flow/gene_flow_analy/script/triplet_mapper.py $speciesTr unbalanced.trp.tsv
#result in unbalanced_triples_raw_count.tre and unbalanced_triples_perc_reticulation_index.tre




############################### 3. calculate the relative contributions of ILS, GTEE, and IH ######
# claculate the gene tree variation
iqtree -t species.tre --gcf BSgenetree.trees --prefix concord

#Final result:
concord.cf.stat


##  GTEE 
#The input sp tree should be the species tree without support values
nw_topology -I rosa_Astral_species_rt.tre > rosa_Astral_species_rt_nosupp.tre

seq-gen -mGTR -l750 -n100  < rosa_Astral_species_rt_nosupp.tre > simulate_tree.fa  


  seq-gen -mGTR -l750 -n100 -a 0.5 -g 4 -z $RANDOM < fixed_tree.tre > simulate_tree.fa  

# separate the simulate seqs into 20 files,becaues we have 20 tips
awk '/^>/{c++} {print > "simulate_" int((c-1)/20)+1 ".fasta"}' simulate_tree.fa

for name in seq/*.fasta;do
  iqtree -s $name -m MFP -B 1000 --bnni -T 50
done

raxmlHPC -f b -t species.tre -z sim_gene.trees -m GTRGAMMA -n ERR

#file result:
RAxML_bipartitions.ERR


# ILS
iqtree -s rosa_tribe_supermatrix.fasta -g rosa_tribe_Astral_species.rt.tre -m MFP -B 1000 --bnni -T 50 -pre rosa_re

#R
library(ape)

# 读取两棵树
tree_mut <- read.tree("rosa_re_rt.treefile")   # IQTREE或RAxML的
tree_coal <- read.tree("rosa_Astral_species_rt.tre")  # ASTRAL或MPEST的

# 提取边和分支长度
edge_mut <- tree_mut$edge
edge_coal <- tree_coal$edge

length_mut <- tree_mut$edge.length
length_coal <- tree_coal$edge.length

# 检查边是否对齐
if (!all(edge_mut == edge_coal)) {
  stop("ERROR: The edge matrices do not match. You must fix the tree order.")
}

# 对齐边：找出每条边（子节点）在mut树中的索引
match_idx <- match(tree_mut$edge[,2], tree_coal$edge[,2])

# 注意：必须对 child（第二列）对齐，因为 tip 和内部节点编号是一致的

# 排序后的 coalescent 分支长度
length_coal_sorted <- tree_coal$edge.length[match_idx]

# θ 计算
theta <- tree_mut$edge.length / length_coal_sorted

# 把每条边的起点、终点和theta输出
df <- data.frame(
  parent = tree_mut$edge[,1],
  child = tree_mut$edge[,2],
  mutation_unit = tree_mut$edge.length,
  coalescent_unit = length_coal_sorted,
  theta = theta
)

# 可选：过滤掉 tip 边，只保留内部节点
Ntip <- length(tree_mut$tip.label)
df_internal <- df[df$child > Ntip, ]

# 写出结果
write.table(df_internal, "theta_per_node.tsv", sep="\t", quote=FALSE, row.names=FALSE)





# IH/Gene flow
library(ape)

# 读取树
tree <- read.tree("unbalancedTriplet.perc.tre")

# 提取基本信息
Ntip <- length(tree$tip.label)
Nnode <- tree$Nnode

# 内部节点的真实编号是：
internal_node_ids <- (Ntip + 1):(Ntip + Nnode)

# 提取对应的 node.label（顺序一一对应）
df <- data.frame(
  node = internal_node_ids,
  reticulate_index = as.numeric(tree$node.label)
)

# 写出为表格
write.table(df, file = "GTEE_err_reticulate_index.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
