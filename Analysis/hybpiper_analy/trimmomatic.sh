#!/bin/bash
# Author: Tao Xiong
# Date: 2025-03.26

input_dir=/home/xiongtao/data/project/Rosaceae/seqdata/final_fastq/part
output_dir=/home/xiongtao/data/project/Rosaceae/seqdata/final_fastq/trimmomatic
#specify the parameter and pathway of the Trimmomatic
trimmomatic_path=/home/xueqin/data/software/Trimmomatic-0.38/trimmomatic-0.38.jar
trimmomatic_params="ILLUMINACLIP:/home/xueqin/data/software/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

#deal with the sample_1.fastq
for sample in ${input_dir}/*_1.fastq; do
#extract the sample name(it means the delete the "_1.fastq")
  sample_name=$(basename ${sample} _1.fastq)
#the left is the sample_1.fastq
  left_reads=${sample}
#the right equal sample_2.fastq
  right_reads=${input_dir}/${sample_name}_2.fastq
  output_prefix=${output_dir}/${sample_name}
#Trimmomatic the data
  java -jar ${trimmomatic_path} PE -threads 15 -phred33 ${left_reads} ${right_reads} ${output_prefix}_1.fq.gz ${output_prefix}_1_unpaired.fastq.gz ${output_prefix}_2.fq.gz ${output_prefix}_2_unpaired.fastq.gz ${trimmomatic_params}
done



#!/bin/bash
input_dir=/home/xiongtao/data/project/Rosaceae/seqdata/final_fastq/single
output_dir=/home/xiongtao/data/project/Rosaceae/seqdata/final_fastq/trimmomatic
#specify the parameter and pathway of the Trimmomatic
trimmomatic_path=/home/xueqin/data/software/Trimmomatic-0.38/trimmomatic-0.38.jar
trimmomatic_params="ILLUMINACLIP:/home/xueqin/data/software/Trimmomatic-0.38/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
for sample in ${input_dir}/*.fastq; do
  sample_name=$(basename ${sample} .fastq)
  reads=${sample}
  output_prefix=${output_dir}/${sample_name}
  java -jar ${trimmomatic_path} SE -threads 15 -phred33 ${reads} ${output_prefix}.fq.gz ${trimmomatic_params}
done