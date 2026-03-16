#!/bin/bash
# Author: Tao Xiong
# Date: 2026-03-15
# Description: 
# ==== 主体代码开始 ====
mkdir -p results

for i in {01..20}; do
    echo "Processing batch $i..."
    
    # 1. 随机抽取 75 个文件路径，存入临时变量
    files=$(ls ../02-trimal/*.fasta | shuf -n 75)
    
    # 2. 使用 pxcat 合并
    # -s 指定输入文件列表，$files 会被展开为多个路径
    # -p 输出分区表, -o 输出超级矩阵
    pxcat -s $files -p results/part_$i.txt -o results/super_$i.fasta
done

parallel 'iqtree -s {} -m MFP -B 1000 --bnni -T 8 -pre {.}' ::: ./*.fasta