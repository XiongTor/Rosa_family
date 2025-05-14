#!/bin/bash
# Author: Tao Xiong
# Date: 2025-05-13
# Description: This script is used to run the tribe ILS flow
# ==== 主体代码开始 ====

##第一次尝试，数据位于：/home/xiongtao/data/project/Rosaceae/ILS/tribe
##时间为2025-05-13

#############################数据准备，hybpiper未筛选直系同源的353数据集合，每个族取了1~3个代表##############################
#mafft
  mkdir mafft trimal
  for name in `ls hybpiper_ILS_result`;do
    echo $(basename $name .FNA).mafft
    mafft --auto hybpiper_ILS_result/$name > mafft/$(basename $name .FNA).mafft
  done

#trimal
  for name in ls mafft/*.mafft; do
    echo ${name}.tri.fasta
    output_file="trimal/$(basename "$name" .mafft).tri.fasta"
    trimal -in "$name" -out "$output_file" -automated1
  done

#build iqtree
for name in trimal/*.tri.fasta; do
  iqtree -s $name -m MFP -B 1000 --bnni -T 15
done

#ASTRAL tree
