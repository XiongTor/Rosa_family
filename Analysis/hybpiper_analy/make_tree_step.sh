# !/usr/bin/bash
# Author: Tao Xiong
# Date: 2025-04.04
# This script is used to record the steps of making tree for hybpiper analysis.

mkdir mafft trimal

for name in *.fasta;do
    echo $(basename $name .fasta).mafft
    mafft --auto $name > mafft/$(basename $name .fasta).mafft
done

for name in ls mafft/*.mafft; do
    tt=$(basename $name .mafft)
    trimal -in $name -out trimal_auto/${tt}.tri.fasta -automated1
done

for name in ls mafft/*.mafft; do
   tt=$(basename $name .mafft)
   trimal -in $name -out trimal/${tt}.tri.fasta -gt 0.8 -st 0.001
done


#prank
python /data/xiongtao/project/Rosaceae/tree/hybpiper_orgin/prank_wrapper.py ./mafft ./prank mafft dna

#genetree
for name in trimal_auto/*.tri.fasta; do
  iqtree -s $name -m MFP -B 1000 --bnni -T 10
done

#connection
pxcat -s trimal_auto/*.fasta -p connection_tree_auto/rosa_sp_partition.txt -o connection_tree_auto/rosa_sp_supermatrix.fasta

SRR27599619
SRR28089552


iqtree -s rosa_sp_supermatrix.fasta --seed 12345 -o SRR27599619,SRR28089552,Elaeagnus_pungens_hangzhou -B 1000 -T 10 -p rosa_sp_partition.txt -m MFP -alrt 1000