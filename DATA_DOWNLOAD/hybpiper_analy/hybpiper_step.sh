# !/usr/bin/bash
# Author: Tao Xiong
# Date: 2025-03.26

#分步进行hybpiper运行
# ```{bash, eval=FALSE, highlight=FALSE}
# **hybpiper assemble**
hybpiper assemble -t_dna Reference_353.fasta -r trimmomatic/*.fq.gz --prefix rosa --bwa --hybpiper_output ./hybpiper
# The parent output directory if supplied using the parameter --hybpiper_output or -o.
#-t_dna目标文件是核苷酸，-t_aa是氨基酸
#-t_aa默认是

while read -r name; do
  hybpiper assemble -t_dna Reference_353.fasta -r ./trimmomatic/$name*.fq.gz --prefix $name --bwa --hybpiper_output ./hybpiper
done <run_hybpiper_list.txt

#namelist.txt：
#EG30
#EG98
#MWL2
#NZ281
#NZ728
#NZ739
#NZ866
#NZ874
#NZ911
#--bwa针对参考序列为核苷酸，BLASTx和DIAMOND针对参考序列为氨基酸，其中BLASTx较慢，DIAMOND可以作为其替代品，--diamond_sensitivity
#**Summary statistics**
ls | grep "_" >namelist.txt
hybpiper stats -t_dna ../Reference_353.fasta gene namelist.txt

#默认输出表格*.tsv（如seq_lengths.tsv），使用--seq_lengths_filename可以进行修改。注：即使参考序列是蛋白文件（amino-acid sequences），该部分还是计算核苷酸（nucleotides）的长度
# **Visualizing results**
hybpiper recovery_heatmap seq_lengths.tsv
#sample_text_size和gene_text_size可以调整文中X、Y轴的文字大小
# **hybpiper retrieve_sequences**
hybpiper retrieve_sequences -t_dna ../Reference_353.fasta --sample_names namelist.txt --fasta_dir ../hybpiper_result dna
#按基因名合并所有物种序列，输出每个基因未比对的序列（unaligned fasta files(one per gene)），使用--fasta_dir输出到指定的文件夹


# **hybpiper paralog_retriever**
hybpiper paralog_retriever namelist.txt -t_dna test_targets.fasta
#组装产生多个contigs时则会存在
#HybPiper检测到多个包含长编码序列的contigs时-默认情况下至少75%的基因序列
#paralog_report.tsv：如果> 1，则为潜在的旁系基因