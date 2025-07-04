#!/bin/bash
# Author: Tao Xiong
# Date: 2025-06-25
# Description: This script is used to analyze the relative contributions of ILS, GTEE, and IH to phylogenetic discordance.
# ==== 主体代码开始 ====

################################ 1. building phyloglenetic tree ###################################
# we need to build a species tree and many bootstrap replicates of gene tree. And all of them should be rooted.
# Get 100 bootstrap gene trees, it means we need 353*100=35300 gene trees.
for seq in `ls ../trimal/*.fasta`; do
  iqtree -s $seq -m MFP --bnni -T 40 -b 100 -pre ./gene_trees/$(basename $seq .fasta)
done

# Astral species tree
cat gene_trees/*.treefile > rosa_sp_genes.tre
java -jar /home/xiongtao/data/software/ASTRAL-master/Astral/astral.5.7.8.jar -i rosa_sp_genes.tre --outgroup Zelkova_schneideriana -o rosa_Astral_species.tre

#Rescale branch lengths
raxmlHPC -f e -t rosa_Astral_species.tre -m GTRGAMMA -s ../rosa_sp_supermatrix.fasta -n rosa_astral_species_br.tre

#Root the tree
pxrr -t rosa_astral_species_br -g Zelkova_schneideriana >rosa_astral_species_br_rt.support

# Delete the support messages in speciestree
Rscript clean_supp.R rosa_Astral_species_rt.tre rosa_Astral_species_rt_clean.tre

#Make sure the species tree do not have the zero branch length
Rscript make_treelength_gt_0.R rosa_Astral_species_rt_clean.tre

## And you should cat all genetrees into one file, and so does the bootstrap trees
cat gene_trees/*.genetrees > rosa_gene_trees.tre
cat gene_trees/*.bootstrapped.tre > rosa_gene_trees_bootstrap.tre


################################ 2. calculate the relative contributions of ILS, GTEE, and IH ######################

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
# !!!! Note, do not use the "./BP.sim.genetrees" as the input, remove the "./", just use "BP.sim.genetrees" is ok. !!!!
cat geneTr_sim/*.sim.genetrees>BP_sim.genetrees.txt
cat BP_sim.genetrees.txt |parallel -j 30 'python ./script/triple_frequency_counter.py {} rosa_Astral_species_rt.tre'

# Find significantly unbalanced triplets and map to species tree
# !!! The GitHub orgin script is wrong, the input file should be the genetrees,not the species tree. You can check the `find_unbalanced_triplets.py`
python /data/xiongtao/tree/ILS/gene_flow/gene_flow_analy/script/find_unbalanced_triplets.py $geneTr

#result in unbalanced.trp.tsv
python /data/xiongtao/tree/ILS/gene_flow/gene_flow_analy/script/triplet_mapper.py $speciesTr unbalanced.trp.tsv
#result in unbalanced_triples_raw_count.tre and unbalanced_triples_perc_reticulation_index.tre




############################### 3. calculate the relative contributions of ILS, GTEE, and IH ######
###  claculate the gene tree variation
iqtree -t species.tre --gcf BSgenetree.trees --prefix concord
#Final result:
concord.cf.stat


###   GTEE 
# use the simulated gene trees to estimate the GTEE error
#The input sp tree should be the species tree without support values
seq-gen -mGTR -l750 -n100 -a 0.5 -g 4 -z $RANDOM < fixed_tree.tre > simulate_tree.fa  

# separate the simulate seqs into 20 files,becaues we have 20 tips
awk '/^>/{c++} {print > "simulate_" int((c-1)/20)+1 ".fasta"}' simulate_tree.fa

for name in seq/*.fasta;do
  iqtree -s $name -m MFP -B 1000 --bnni -T 50
done

cat seq/*.treefile > sim_gene.trees
raxmlHPC -f b -t species.tre -z sim_gene.trees -m GTRGAMMA -n ERR

#file result: RAxML_bipartitions.ERR
Rscript ~/data/scripts/ILS_GTEE_IH/get_node_inf.R RAxML_bipartitions.ERR

###  ILS
iqtree -s rosa_tribe_supermatrix.fasta -g rosa_tribe_Astral_species.rt.tre -m MFP -B 1000 --bnni -T 50 -pre rosa_re

Rscript Rscript ~/data/scripts/ILS_GTEE_IH/get_theta_value.R rosa_Astral_species_rt.tre rosa_re.treefile
#result: theta_per_node.tsv

### IH/Gene flow
Rscript ~/data/scripts/ILS_GTEE_IH/get_node_inf.R unbalanced_triples_perc_reticulation_index.tre

#