#!/bin/bash
# Author: Tao Xiong
# Date: 2026-02-09
# Description: 在完成trimal,treeshrink等一系列步骤后，统计最终序列文件中的每个物种在多少个基因序列文件中仍有数据
# ==== 主体代码开始 ====
# 1. 提取所有文件的 ID 并计数
seqkit seq --name *.fasta | cut -d " " -f 1 | sort | uniq -c | sort -nr | sed 's/^[[:blank:]]*//;s/[[:blank:]]\+/,/' > species_gene_counts.csv

# 2. 赋予列名
sed -i '1i Count,Species' species_gene_counts.csv

# --------------------------------------------------------------------------------------------
# 在完成R比对，得到族的名字后，开始下列步骤

sed -i 1d final_to_orthofinder_built_tree.csv
# 1. 创建存放结果的文件夹
mkdir -p selected_species_genes

# 2. 批量提取循环
# 假设你的原始序列在当前目录下，后缀为 .fasta
for file in *.fasta; do
    # 使用 seqkit grep 根据名单匹配 ID
    # -f 指定名单文件
    # -o 指定输出文件名
    seqkit grep -f final_to_orthofinder_built_tree.csv "$file" -o "selected_species_genes/${file%.fasta}_extracted.fasta"
done

echo "提取完成！结果保存在 selected_species_genes 目录下。"


mv selected_species_genes ./seq_data
mkdir genetrees

for name in seq_data/*.fasta;do
  tt=$(basename $name .fasta)
  iqtree -s $name -m MFP -B 1000 --bnni -T 30 -pre ./genetrees/${tt}
done