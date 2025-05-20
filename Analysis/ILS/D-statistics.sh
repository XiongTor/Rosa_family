#!/bin/bash
# Author: Tao Xiong
# Date: 2025-05-16
# Description: This script is used to calculate the D-statistics to check the gene introgression.
# ==== 主体代码开始 ====

########################## get the SNPs VCF file ################################
# install snp-site
# you can see the website of GitHub:https://github.com/sanger-pathogens/snp-sites
   apt-get install snp-sites

# prepare the input files
# the input files is the fasta file you used to build the phylogenetic tree, it shoule be a supermatrix file.
pxcat -s ../trimal/*.fasta -p rosa_tribe_partition.txt -o rosa_tribe_supermatrix.fasta

# get the SNPs VCF file
snp-sites -mvp -o rosa_snp rosa_sp_supermatrix.fasta

########################## D-statistics calculation ##############################
Dsuite Dtrios rosa_snp.vcf sets.txt -t rosa_Astral_species.rt.nolength.tre -o sample

# 计算f-branch值
Dsuite Fbranch rosa_Astral_species.rt.nolength.tre sets_tree.txt >fbranch.out

# 用dtools.py脚本绘制f-branch图
~/data/software/Dsuite/utils/dtools.py fbranch.out rosa_Astral_species.rt.nolength.tre --outgroup Elaeagnus_pungens_hangzhou --use_distances --dpi 400 --tree-label-size 15
