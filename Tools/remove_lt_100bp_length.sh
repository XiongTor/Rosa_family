#!/bin/bash
# Author: Tao Xiong
# Date: 2025-10-11
# Description: used to remove the sequence length less than 100bp
# ==== 主体代码开始 ====
#寻找序列长度短于100bp的序列并移除
n=0
for name in ./RLWP/*.fasta;do
  n=$((n+1))
  echo $n
  tt=$(seqkit fx2tab -l -n $name|cut -f2|datamash median 1)
  if [[ $tt -le 100 ]];then
    echo $name
    echo $name>>RLWP_less100.txt
    mv $name ./less100_length/RLWP
  fi
done


#剔除序列长度短于100bp的基因序列构建的基因树
while read -r line;do
  tt=$(echo $line|cut -f3 -d'/'|cut -f1 -d'.')
  mv 07-genetrees/RT/${tt}.* 08-genetree_less100/RT
done<RT_less100.txt

  tt=$(cat $line|cut -f3 -d'/'|cut -f1 -d'.')
