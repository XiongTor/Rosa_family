#!/bin/bash
# Author: Tao Xiong
# Date: 2025-10-11
# Description: used to remove the sequence length less than 100bp
# ==== 主体代码开始 ====
#寻找序列长度短于100bp的序列并移除

n=0
for name in ./trimal/*.fasta;do
  n=$((n+1))
  echo $n
  tt=$(seqkit fx2tab -l -n $name|cut -f2|datamash median 1)
  if [[ $tt -le 100 ]];then
    echo $name
    echo $name>>trimal_less100.txt
    mv $name ./remove_less100_trimal
  fi
done


#剔除序列长度短于100bp的基因序列构建的基因树
while read -r line;do
  tt=$(echo $line|cut -f3 -d'/'|cut -f1,2 -d'_')
  mv genetrees_python/${tt}* less100_after_treeshrink_genetree_python
done<trimal_less100.txt


#剔除哪些序列中全为gap或者N的序列
while read -r name;do
  echo $name
  seqkit grep -s -r -p "[ACGTacgt]" rm_toomuch_gap_trimal/${name}.tri.fasta> 02-trimal/${name}_trimmed_haveoutg_output.tri.fasta
done<rebulid_list.txt