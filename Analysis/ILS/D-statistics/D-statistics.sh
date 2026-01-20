#!/bin/bash
# Author: Tao Xiong
# Date: 2025-05-16
# Description: This script is used to calculate the D-statistics to check the gene introgression.

# ==== 主体代码开始 ====

#chatgpt获得vcf文件
cov=(0.04,0.5,0.6,0.7,0.8,0.9)
maf=(0.01,0.05)



python ~/data/scripts/msa2vcf_biallelic.py ../rosa_MO_orthofinder_supermatrix.fasta rosa_orthofinder.vcf --min_cov 0.8 --min_maf 0.05 

########################## D-statistics calculation ##############################

Dsuite Dtrios rosa_orthofinder.vcf ../sets.txt -t ../rosa_orthofinder_MO_treeshrink_sp_rt.tre -o result

ruby ~/data/scripts/Dsuite/plot_d.rb result_tree.txt ../plot_order.txt 0.7 species_sets_no_geneflow_BBAA_D.svg

ruby ~/data/scripts/Dsuite/plot_f4ratio.rb result_tree.txt ../plot_order.txt 0.2 species_sets_no_geneflow_BBAA_f4ratio.svg  

# 计算f-branch值

Dsuite Fbranch ../rosa_orthofinder_MO_treeshrink_sp_rt.tre result_tree.txt >fbranch.out
  

# 用dtools.py脚本绘制f-branch图
dtools.py fbranch.out ../rosa_orthofinder_MO_treeshrink_sp_rt.tre --outgroup Zelkova_schneideriana --use_distances --dpi 500 --tree-label-size 20 