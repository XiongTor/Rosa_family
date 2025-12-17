#!/bin/bash
# Author: Tao Xiong
# Date: 2025-12-16
# Description: This script 
# ==== 主体代码开始 ====

mkdir -p tmp

# 获取迭代次数（列数）
for iter in $(seq 1 20); do
    echo "Processing Iteration ${iter}..."
    gene_list=$(awk -F',' -v col=$((iter)) 'NR>1 {gsub(/\r/, "", $col); print $col}' treeness_Chloroplast_like_nuclear_gene_genelist.csv)
    for name in ${gene_list}; do
        cp ../02-trimal/${name}_* tmp
    done
   pxcat -s tmp/*.fasta -p Iteration_${iter}_partition.txt -o Iteration_${iter}_supermatrix.fasta
   rm tmp/*
done

mkdir Chloroplast_like
mv *.fasta *.txt Chloroplast_like


iqtree -s Chloroplast_like/$fasta -m MFP -B 1000 --bnni -T 5 -pre ./Chloroplast_like

mkdir treefile
cp Chloroplast_like/*.treefile treefile

for name in ./*.treefile;do
  mm=$(basename $name .treefile)
  echo $mm
  pxrr -t $name -g Elaeagnus_angustifolia,Zelkova_schneideriana,Morus_indica > ./${mm}_rt.tre
  pxrlt -t ${mm}_rt.tre -c ../../old_name.txt -n ../../new_name.txt>./${mm}_rt_rn.tre
  Rscript ~/data/scripts/mono.R ./${mm}_rt_rn.tre
done

#!/bin/bash

parallel 'iqtree -s {} -m MFP -B 1000 --bnni -T 5 -pre {.}' ::: Mean_evolutionary_rate/Chloroplast_like/*.fasta
parallel 'iqtree -s {} -m MFP -B 1000 --bnni -T 5 -pre {.}' ::: Proportion_parsimony_informative/Chloroplast_like/*.fasta
parallel 'iqtree -s {} -m MFP -B 1000 --bnni -T 5 -pre {.}' ::: Proportion_variable_sites/Chloroplast_like/*.fasta