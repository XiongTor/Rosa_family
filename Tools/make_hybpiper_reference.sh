# !/usr/bin/bash
# Author: Tao Xiong
# Date: 2025-04.03
# Description: This script is used to get the reference sequence

#1.首先将easy353下载的参考序列进行合并（/data/xueqin/Project/Prunus/Prunus_353_easy353/353_ref_Prunus/353gene）
cat *.fasta >../../con_353.fasta
#2.命名为物种-基因名
sed -i 's/>\([0-9]*\) Gene_Name:[^ ]* Species:\([^ ]*\) Repository:[^ ]* Sequence_ID:[^ ]*/>\2-\1/g' con_353.fasta
#3.删除合并序列中的重复信息（>物种名称+序列两行信息）
awk 'BEGIN { RS=">"; FS="\n" } NR>1 && seen[$1]++ { next } { printf ">%s", $0 }' con_353.fasta > > con_norepeat_353_file.fa #修改参考序列格式（按照每个物种排列，且物种与基因之间空格，且60个字符换行，如）
sed -E -i 's/>>>/>/g'  con_norepeat_353_file.fa #删去出现>>的情况
#3.# 将序列转换成一行 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' con_353.fasta | sed 1d - > con_norepeat_353_file.fa
#4.修改格式为无重复，每60个字符换行
#物种名提取并排序
grep ">" con_norepeat_353_file.fa | sort | uniq > name2.txt
# 提取物种名及相关序列
while read -r line; do  grep -A 1 "$line" con_norepeat_353_file.fa >> out.txt ; done < name2.txt
# 60个字符分隔
awk '!/^>/ { gsub(/.{60}/,"&\n"); print; next } { print }' out.txt > out1.txt
#去除空行信息 
grep -v '^$' out1.txt > con_353_finall.fa