#! /bin/bash
#本脚本用于合并hybpiper和paftol的353数据，基本原则是将paftol的数据经过处理后插入到hybpiper中
#同时为保证格式一致，最终结果的title中均只保留物种名

#指定hybpiper和paftol所在的文件夹
hybpiper=$1
paftol=$2

## 针对hybpiper进行处理 ##===========

#创建新的文件夹
mkdir comb_result_genes

#去掉hybpiper文件中的换行符
for name in $hybpiper/*.FNA;do
    tt=$(basename $name|cut -f1 -d'_')
    awk '/^>/ {printf("\n%s\n",$0);next;} { printf("%s",$0);} END {printf("\n");}' $name |sed '1d'> comb_result_genes/${tt}.fasta
 done

## 针对paftol进行处理 ##===========

#去掉paftol文件中的换行符,提取对应序列并放入到comb_result_genes对应的基因序列文件中
for name in $paftol/*;do
    grep ">" $name | awk -F'[ :]' '{print $1}'|sed 's/>//g'>./${name}_gene_list.txt      #提取一个物种的全部353基因名称
    while read -r gene;do
        #去除换行符
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' $name| \
        #提取基因名和物种名
        sed -E 's/>([0-9]+) Gene_Name:[^ ]+ Species:([^ ]+) Repository:[^ ]+ Sequence_ID:[^ ]+/>\2 single_hit \1/g'| \
        #抓取到对应基因序列并插入到comb_result_genes对应的文件中,同时去除基因名,加上single_hit
        grep -A 1 -w $gene|awk -F' ' '{print $1, $2}' >> comb_result_genes/${gene}.fasta
    done<./${name}_gene_list.txt
    rm ./${name}_gene_list.txt
done

##将所有的序列转化为60碱基一换行的格式
mkdir comb_result_genes_60
for seq in comb_result_genes/*.fasta;do
    tt=$(basename $seq)
    seqkit seq -w 60 $seq > comb_result_genes_60/$tt
done