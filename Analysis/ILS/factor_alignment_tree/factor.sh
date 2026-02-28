#!/bin/bash
# Author: Tao Xiong
# Date: 2025-12-07
# Description: Use phykit and AMAS to get alignment and tree statistics
# ==== 主体代码开始 ====
trimal=$1
genetree_reroot=$2

#####################   use AMAS to get alignment statistics #####################
AMAS.py summary -f fasta -d dna -i $trimal/*.fasta -c 12

#get Alignment_name	No_of_taxa	Alignment_length	Total_matrix_cells	Missing_percent	No_variable_sites	Proportion_variable_sites	Proportion_parsimony_informative	GC_content	N	-	?
cut -f1,2,3,4,6,7,8,10,12,28,30,31 summary.txt>summary_final.txt

echo "finish getting alignment statistics"
#####################   use phykit   #################################################
# get evolutionary rate per site
i=0
for file in $trimal/*.fasta; do
    i=$((i+1))
    echo $i
    name=$(basename $file)
	tt=$(phykit evo_rate_per_site $file|datamash mean 2)
	tt_without0=$(phykit evo_rate_per_site $file| awk '$2>0'| datamash mean 2)
	echo -e "${name}\t${tt}\t${tt_without0}" >> evo_rate.txt
done

sed -i '1i Alignment_name\tMean_evolutionary_rate\tMean_evolutionary_rate_without0' evo_rate.txt

echo "1.finish getting evolutionary rate per site"

# average bipartition support
i=0
for file in $genetree_reroot/*.tre; do
    i=$((i+1))
    echo $i
    name=$(basename $file)
	tt=$(phykit bipartition_support_stats $file -v|cut -f1 -d' '|datamash mean 1)
	echo -e "${name}\t${tt}" >> average_bipartition.txt
done
sed -i '1i Alignment_name\tMean_bipartition' average_bipartition.txt

echo "2.finish getting average bipartition support"

# average teminal branch length
i=0
for file in $genetree_reroot/*.tre; do
    i=$((i+1))
    echo $i
    name=$(basename $file)
	tt=$(phykit terminal_branch_stats $file -v|cut -f1 -d' '|datamash mean 1)
	echo -e "${name}\t${tt}" >> teminal_branch_length.txt
done
sed -i '1i Alignment_name\tMean_teminal_branch' teminal_branch_length.txt

echo "3.finish getting average teminal branch length"

# average internal branch length
i=0
for file in $genetree_reroot/*.tre; do
    i=$((i+1))
    echo $i
    name=$(basename $file)
	tt=$(phykit internal_branch_stats $file -v|cut -f1 -d' '|datamash mean 1)
	echo -e "${name}\t${tt}" >> internal_branch_length.txt
done
sed -i '1i Alignment_name\tMean_internal_branch' internal_branch_length.txt

echo "4.finish getting average internal branch length"

# tree_ness
i=0
for file in $genetree_reroot/*.tre; do
    i=$((i+1))
    echo $i
    name=$(basename $file)
	tt=$(phykit treeness $file)
	echo -e "${name}\t${tt}" >> treeness.txt
done
sed -i '1i Alignment_name\tTreeness' treeness.txt

echo "5.finish getting tree_ness"
# amino acid substitution saturation estimated

i=0
for file in $genetree_reroot/*.tre; do
    i=$((i+1))
    echo $i
    name=$(basename $file _rt.tre)
	tt=$(phykit sat -a $trimal/${name}.fasta -t $file)
	echo -e "${name}\t${tt}" >> saturation.txt
done
sed -i '1i Alignment_name\tsaturation\tabsolute_saturation' saturation.txt

echo "6.finish getting amino acid substitution saturation estimated"

# Relative composition variability
i=0
for file in $trimal/*.fasta; do
    i=$((i+1))
    echo $i
	name=$(basename $file)
	tt=$(phykit relative_composition_variability $file)
	echo -e "${name}\t${tt}" >> RCV.txt
done
sed -i '1i Alignment_name\tRCV' RCV.txt

echo "7.finish getting Relative composition variability"

# Treeness over RCV
i=0
for file in $genetree_reroot/*.tre; do
    i=$((i+1))
    echo $i
    name=$(basename $file .rt.tre)
	tt=$(phykit treeness_over_rcv -a $trimal/${name}.fasta -t $file)
	echo -e "${name}\t${tt}" >> Treeness_over_RCV.txt
done
sed -i '1i Alignment_name\ttreeness/RCV\ttreeness\tRCV' Treeness_over_RCV.txt
echo "8.finish getting Treeness over RCV"



i=0
while read -r file;do
    i=$((i+1))
    echo $i
    name=$(basename $file _mad.tre)
	tt=$(phykit treeness_over_rcv -a 02-trimal/${name}.fasta -t $file)
	echo -e "${name}\t${tt}" >> Treeness_over_RCV.txt
done<mad_tree_list_for_Treeness.txt