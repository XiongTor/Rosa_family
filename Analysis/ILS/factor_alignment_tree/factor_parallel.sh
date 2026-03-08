#!/bin/bash
trimal=$1
genetree_reroot=$2

#####################  AMAS alignment statistics  #####################
AMAS.py summary -f fasta -d dna -i $trimal/*.fasta -c 12
cut -f1,2,3,4,6,7,8,10,12,28,30,31 summary.txt > summary_final.txt
echo "finish getting alignment statistics"

#####################  并行运行所有循环块  #####################

# 1. evolutionary rate per site
{
    for file in $trimal/*.fasta; do
        name=$(basename $file)
        tt=$(phykit evo_rate_per_site $file | datamash mean 2)
        tt_without0=$(phykit evo_rate_per_site $file | awk '$2>0' | datamash mean 2)
        echo -e "${name}\t${tt}\t${tt_without0}" >> evo_rate.txt
    done
    sed -i '1i Alignment_name\tMean_evolutionary_rate\tMean_evolutionary_rate_without0' evo_rate.txt
    echo "1. finish getting evolutionary rate per site"
} &

# 2. average bipartition support
{
    for file in $genetree_reroot/*.tre; do
        name=$(basename $file)
        tt=$(phykit bipartition_support_stats $file -v | cut -f1 -d' ' | datamash mean 1)
        echo -e "${name}\t${tt}" >> average_bipartition.txt
    done
    sed -i '1i Alignment_name\tMean_bipartition' average_bipartition.txt
    echo "2. finish getting average bipartition support"
} &

# 3. terminal branch length
{
    for file in $genetree_reroot/*.tre; do
        name=$(basename $file)
        tt=$(phykit terminal_branch_stats $file -v | cut -f1 -d' ' | datamash mean 1)
        echo -e "${name}\t${tt}" >> teminal_branch_length.txt
    done
    sed -i '1i Alignment_name\tMean_teminal_branch' teminal_branch_length.txt
    echo "3. finish getting average terminal branch length"
} &

# 4. internal branch length
{
    for file in $genetree_reroot/*.tre; do
        name=$(basename $file)
        tt=$(phykit internal_branch_stats $file -v | cut -f1 -d' ' | datamash mean 1)
        echo -e "${name}\t${tt}" >> internal_branch_length.txt
    done
    sed -i '1i Alignment_name\tMean_internal_branch' internal_branch_length.txt
    echo "4. finish getting average internal branch length"
} &

# 5. treeness
{
    for file in $genetree_reroot/*.tre; do
        name=$(basename $file)
        tt=$(phykit treeness $file)
        echo -e "${name}\t${tt}" >> treeness.txt
    done
    sed -i '1i Alignment_name\tTreeness' treeness.txt
    echo "5. finish getting treeness"
} &

# 6. saturation
{
    for file in $genetree_reroot/*.tre; do
        name=$(basename $file .rt.tre)
        tt=$(phykit sat -a $trimal/${name}.fasta -t $file)
        echo -e "${name}\t${tt}" >> saturation.txt
    done
    sed -i '1i Alignment_name\tsaturation\tabsolute_saturation' saturation.txt
    echo "6. finish getting saturation"
} &

# 7. RCV
{
    for file in $trimal/*.fasta; do
        name=$(basename $file)
        tt=$(phykit relative_composition_variability $file)
        echo -e "${name}\t${tt}" >> RCV.txt
    done
    sed -i '1i Alignment_name\tRCV' RCV.txt
    echo "7. finish getting Relative composition variability"
} &

# 8. Treeness over RCV
{
    for file in $genetree_reroot/*.tre; do
        name=$(basename $file .rt.tre)
        tt=$(phykit treeness_over_rcv -a $trimal/${name}.fasta -t $file)
        echo -e "${name}\t${tt}" >> Treeness_over_RCV.txt
    done
    sed -i '1i Alignment_name\ttreeness/RCV\ttreeness\tRCV' Treeness_over_RCV.txt
    echo "8. finish getting Treeness over RCV"
} &

# 等待所有后台任务完成
wait
echo "All tasks finished!"