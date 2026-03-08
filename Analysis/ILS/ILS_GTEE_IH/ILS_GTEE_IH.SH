#!/bin/bash
# Author: Tao Xiong
# Date: 2025-06-25
# Description: This script is used to analyze the relative contributions of ILS, GTEE, and IH to phylogenetic discordance.
# ==== 主体代码开始 ====

################################ 0. prepare the input files ###################################
Input:

$speciesTr：One rooted best-estimated species tree with branch lengths measured in coalescent units (can be estimated from MPEST/ASTRAL on ML gene trees, newick format);
$speciesTr_BP：Rooted bootstrap species trees with branch lengths measured in coalescent units (can be obtained by running MPEST/ASTRAL on bootstrap gene trees, newick format);
$geneTr：Rooted gene trees, newick format.

################################ 1. building phyloglenetic tree ###################################
# we need to build a species tree and many bootstrap replicates of gene tree. And all of them should be rooted.
# Get 100 bootstrap gene trees, it means we need 353*100=35300 gene trees.

# For species tree
for seq in `ls ../trimal/*.fasta`; do
  iqtree -s $seq -m MFP --bnni -T 40 -b 100 -pre ./gene_trees/$(basename $seq .fasta)
done

# Astral species tree
cat gene_trees/*.treefile > rosa_sp_genes.tre
java -jar /home/xiongtao/data/software/ASTRAL-master/Astral/astral.5.7.8.jar -i rosa_sp_genes.tre --outgroup Zelkova_schneideriana -o rosa_Astral_species.tre

#Rescale branch lengths,for count the theta value
raxmlHPC -f e -t rosa_Astral_species.tre -m GTRGAMMA -s ../rosa_sp_supermatrix.fasta -n rosa_astral_species_br.tre

#Root the tree
pxrr -t rosa_astral_species_br -g Zelkova_schneideriana >rosa_astral_species_br_rt.support

# Delete the support messages in speciestree
Rscript clean_supp.R rosa_Astral_species_rt.tre rosa_Astral_species_rt_clean.tre

#Make sure the species tree do not have the zero branch length
Rscript make_treelength_gt_0.R rosa_Astral_species_rt_clean.tre

# For bootstrap gene trees
#若IQTEE中采用的是-B，采用了近似随机1000次采样的bootstrap值计算，则随机采样100棵树并合并bootstrap树的集合
for file in ../03-genetrees/*.ufboot; do
  tt=$(basename $file .ufboot)
  shuf -n 100 "$file" > ./bootstrap_100/${tt}_100.ufboot
done

#对每个基因的100个bootstrap树集合进行置根
for tree in ./*.ufboot; do
    tt=$(nw_labels -I $tree|sort|uniq|grep -E "Elaeagnus_angustifolia|Zelkova_schneideriana|Morus_indica"|paste -sd,)
    mm=$(basename $tree .ufboot)
    pxrr -t $tree -g $tt > ./${mm}_rt.tre
done

# 合并每个基因的bootstrap树的第i行，生成第i个bootstrap replicate的基因树集合
# 每个文件有 100 行
for i in $(seq 1 100); do
    out="bootstrap_${i}.tre"
    echo "生成 $out"
    for f in ../bootstrap_100/*.ufboot; do
        sed -n "${i}p" "$f"
    done > "$out"
done

#Build species tree for each bootstrap replicate by ASTRAL-Pro3
for i in $(seq 1 100); do
  echo "======== start bootstral_${i} ========"
  astral-pro3 -i bootstrap_${i}.tre -o rosa_ags353_bootstrap_${i}_sp.tre
  pxrr -t rosa_ags353_bootstrap_${i}_sp.tre -g Elaeagnus_angustifolia,Zelkova_schneideriana,Morus_indica>rosa_ags353_bootstrap_${i}_sp_rt.tre
  Rscript ../script/clean_supp.R rosa_ags353_bootstrap_${i}_sp_rt.tre rosa_ags353_bootstrap_${i}_sp_rt_nosupp.tre
done

## And you should cat all genetrees into one file, and so does the bootstrap trees
cat gene_trees/*.genetrees > rosa_gene_trees.tre

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
iqtree -t rosa_ags353_treeshrink_sp_rt_nosupp_nozero_oneoutg_final.tre --gcf BSgenetree.trees --prefix concord
#BSgenetree.trees is the bootstrap gene trees
#Final result:
concord.cf.stat


###   GTEE 
# use the simulated gene trees to estimate the GTEE error
#The input sp tree should be the species tree without support values
seq-gen -mGTR -l500 -n317 -a 0.5 -g 4 -z $RANDOM < rosa_ags353_treeshrink_sp_rt_nosupp_nozero_oneoutg_final.tre > simulate_tree.fa  

# separate the simulate seqs into 20 files,becaues we have 20 tips
# you can use the mean /median length of your real seqs to set the length of simulate seqs
awk '/^>/{c++} {print > "simulate_" int((c-1)/80)+1 ".fasta"}' simulate_tree.fa

ls *.fasta | parallel -j 3 "iqtree3 -s {} -m MFP -B 1000 --bnni -T 15"

for name in seq/*.fasta;do
  iqtree -s $name -m MFP -B 1000 --bnni -T 50
done

cat seq/*.treefile > sim_gene.trees
raxmlHPC -f b -t species.tre -z sim_gene.trees -m GTRGAMMA -n ERR

#file result: RAxML_bipartitions.ERR
Rscript ~/data/scripts/ILS_GTEE_IH/get_node_inf.R RAxML_bipartitions.ERR

###  ILS
iqtree -s rosa_tribe_supermatrix.fasta -g rosa_tribe_Astral_species.rt.tre -m MFP -B 1000 --bnni -T 50 -pre rosa_re

Rscript ~/data/scripts/ILS_GTEE_IH/get_theta_value.R rosa_Astral_species_rt.tre rosa_re.treefile
#result: theta_per_node.tsv

### IH/Gene flow
Rscript ~/data/scripts/ILS_GTEE_IH/get_node_inf.R unbalanced_triples_perc_reticulation_index.tre

###  combine the results
# New you should conbine all the results into one file , like this:
# nodeID	genetree_var	ILS		GTEE	IH
# 1			963				0.5		0.5		1
# 2			4				0.5		0.5		11
# 3			801				4		27		4
# 4			849				4		62		100

Rscript /data/xiongtao/scripts/ILS_GTEE_IH/relaimpo.R ILS_GTEE_IH.csv