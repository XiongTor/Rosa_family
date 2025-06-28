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

## simulated the sequnce to evaluate the GTEE LEVEL

#The input sp tree should be the species tree without support values
nw_topology -I rosa_Astral_species_rt.tre > rosa_Astral_species_rt_nosupp.tre

seq-gen -mGTR -l750 -n100  < rosa_Astral_species_rt_nosupp.tre > simulate_tree.fa  

for i in {1..100}; do
  seq-gen -mGTR -l750 -n1 0a-z $RANDOM < rosa_Astral_species_rt_nosupp.tre > simulate_tree_${i}.fa  
done
# separate the simulate seqs into 20 files,becaues we have 20 tips
awk '/^>/{c++} {print > "simulate_" int((c-1)/20)+1 ".fasta"}' simulate_tree.fa

for name in seq/*.fasta;do
  iqtree -s $name -m MFP -B 1000 --bnni -T 15
done