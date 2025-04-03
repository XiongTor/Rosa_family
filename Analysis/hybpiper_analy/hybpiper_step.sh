# !/usr/bin/bash
# Author: Tao Xiong
# Date: 2025-03.26

#分步进行hybpiper运行

# hybpiper assemble
hybpiper assemble -t_dna hybpiper/Reference_353.fasta -r ./final_fastq/trimmomatic/SRR26662633*.fq.gz --prefix Neillia --bwa --hybpiper_output ./Neillia_chloro
# The parent output directory if supplied using the parameter --hybpiper_output or -o.
#-t_dna目标文件是核苷酸，-t_aa是氨基酸
#-t_aa默认是
#--bwa针对参考序列为核苷酸，BLASTx和DIAMOND针对参考序列为氨基酸，其中BLASTx较慢，DIAMOND可以作为其替代品，--diamond_sensitivity

# loop
while read -r name; do
  hybpiper assemble -t_dna Reference_353.fasta -r ./$name*.fq.gz --prefix $name --bwa --hybpiper_output ./hybpiper
done <run_hybpiper_list.txt

# Summary statistics
ls | grep "_" >namelist.txt
hybpiper stats -t_dna ../Reference_353.fasta gene namelist.txt
#默认输出表格*.tsv（如seq_lengths.tsv），使用--seq_lengths_filename可以进行修改。注：即使参考序列是蛋白文件（amino-acid sequences），该部分还是计算核苷酸（nucleotides）的长度

# Visualizing results
hybpiper recovery_heatmap seq_lengths.tsv
#sample_text_size和gene_text_size可以调整文中X、Y轴的文字大小

# hybpiper retrieve_sequences
hybpiper retrieve_sequences -t_dna ../Reference_353.fasta --sample_names namelist.txt --fasta_dir ../hybpiper_result dna
#按基因名合并所有物种序列，输出每个基因未比对的序列（unaligned fasta files(one per gene)），使用--fasta_dir输出到指定的文件夹
