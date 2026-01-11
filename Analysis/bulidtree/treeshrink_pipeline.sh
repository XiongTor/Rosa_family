#!/bin/bash
# Author: Tao Xiong
# Date: 2026-01-11
# Description: This script is to filter gene sequences with length less than 100bp and make treeshrink analysis.
# ==== 主体代码开始 ====



# 1. 过滤序列长度，建议直接查看最小的几个文件进行处理，代码执行有一定问题，速度却不一定更快
# 寻找序列长度短于100bp的序列并移除
n=0
for name in ./00-MO_orgin_seq/*.fasta;do
  n=$((n+1))
  echo $n
  tt=$(seqkit fx2tab -l -n $name|cut -f2|datamash median 1)
  if [[ $tt -le 100 ]];then
    echo $name
    echo $name>>MO_less100.txt
    mv $name ./less100_seq
  fi
done


#剔除序列长度短于100bp的基因序列构建的基因树
while read -r line;do
  tt=$(echo $line|cut -f3 -d'/'|cut -f1 -d'.')
  mv 01-MO_orgin_genetrees/${tt}.* genetree_less100
done<MO_less100.txt

# 2. 进行treeshrink分析

#先将每个物种的序列换行符去掉
n=0
for name in 00-MO_orgin_seq/*.fasta;do
  n=$((n+1))
  echo $n
  tt=$(basename $name .fasta)
  seqkit seq ${name} -w 0 > ./MO_orgin_seq_oneline/${tt}_oneline.fasta
done

# 将外类群序列单独提取出来放到对应文件夹中
m=0
for name in ./MO_orgin_seq_oneline/*.fasta;do
  m=$((m+1))
  echo $m
  tt=$(basename $name _oneline.fasta)
  grep -E -A 1 --no-group-separator "Elaeagnus_angustifolia|Zelkova_schneideriana|Morus_indica" $name> ./MO_outg_seq/${tt}_outgroup.fasta
done


#去除所有基因树中的外类群，并将对应的序列放到对应文件夹中，序列更名为input.fasta，基因树更名为input.tre
m=0
for tree in 01-MO_orgin_genetrees/*.treefile;do
  m=$((m+1))
  echo $m
  tt=$(basename $tree .treefile|cut -f1 -d'.')
  mkdir ./MO_treeshrink/$tt
  pxrmt -t $tree -n Elaeagnus_angustifolia,Zelkova_schneideriana,Morus_indica >./MO_treeshrink/$tt/input.tre
  pxrms -s ./MO_orgin_seq_oneline/${tt}*.fasta -n Elaeagnus_angustifolia,Zelkova_schneideriana,Morus_indica >./MO_treeshrink/$tt/input.fasta
done

#运行treeshrink,注意路径
nohup run_treeshrink.py -i ./ -t input.tre -a input.fasta -b 20 -q 0.2 > ./input.tree.treeshrinklog.txt &

#查看全部的长枝报表，本步骤中阈值参数设置无用，最终输出的报表一致，为所有长枝的相关信息
# run_treeshrink.py  -t MO_orthofinder_genetrees.tre -q 0.01 -b 1 -o treeshrink_multi_MO -O shrunk_0.01

#更改输出文件名称,增加外类群序列
a=0
for name in `ls -d */`;do
  a=$((a+1))
  echo $a
  tt=$(basename $name /)
  mv $tt/output.tre $tt/${tt}_noout_output.tre
  cat ../MO_outg_seq/${tt}*.fasta >>$tt/output.fasta
  mv $tt/output.fasta $tt/${tt}_noout_output.fasta
done


#建树
ls -d */|sed 's/\///g'>namelist.txt


while read -r line;do
  cp ${line}/*output.fasta seq_data 
done<namelist.txt

### 值得注意的是，treeshrink似乎会在去除序列后自动做序列修剪，导致按照此流程，插入的原序列的外类群会和treeshrink后的序列无法对齐，从而导致建树错误，因此建议在此时重新比对修剪


for name in `ls seq_data`;do
  echo $(basename $name .fasta).mafft
  mafft --auto seq_data/$name > mafft/$(basename $name .fasta).mafft
done

for name in ls mafft/*.mafft; do
  echo ${name}.tri.fasta
  output_file="trimal/$(basename "$name" .mafft).tri.fasta"
  trimal -in "$name" -out "$output_file" -automated1
done

for name in trimal/*.fasta;do
  tt=$(basename $name .fasta)
  iqtree -s $name -m MFP -B 1000 --bnni -T 10 -pre ./genetrees/${tt}
done

-----------------
while read -r line;do
  iqtree -s $line/${line}_trimmer_haveout_output.fasta -m MFP -B 1000 --bnni -T 10 -pre ./genetree_b20/${line}
done<namelist.txt