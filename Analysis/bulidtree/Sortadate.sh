# !/usr/bin/bash
# Author: Xiongtao
# 代码位置：/data/xiongtao/scripts/SortaDate/

#Sortadate流程适用于溯祖法，本脚本主要是针对溯祖法进行

# $1 包含有所有置根后的基因树的文件夹
# $2 置根后的物种树
# $3 外类群的名字
# $4 想要筛选出的基因树的个数

#Get the root-to-tip variance with:
python /data/xiongtao/scripts/SortaDate/src/get_var_length.py $1 --flend .tre --outf var --outg $3
# Tree_reroot：为包含所有置根后的基因树的文件夹(注意是最优树)
# Barbeya_oleoides：为使用的外类群的名字

#Get the bipartition support with:
python /data/xiongtao/scripts/SortaDate/src/get_bp_genetrees.py $1 $2 --flend .tre --outf bp
# rosa_Astral_species.tre：为置根后的物种树

#Combine the results from these two runs with 
python /data/xiongtao/scripts/SortaDate/src/combine_results.py var bp --outf comb

#Sort and get the list of the good genes with 
python /data/xiongtao/scripts/SortaDate/src/get_good_genes_xt.py comb --max $4 --outf result.txt

## 此处使用的get_good_genes_xt.py是重新编写的脚本，主要是可以指定多个外类群

# 提取出筛选出的基因树的名字
cat result.txt|sed '1d'|cut -f1 -d' '>tree_name.txt

# 提取出对应的基因的名字
cat tree_name.txt|cut -f2 -d'.'>gene_sortadate_name.txt

mkdir tree_sortadate

# 将对应的基因树复制到新文件夹中
parallel cp $1/{} tree_sortadate/{} <tree_name.txt











