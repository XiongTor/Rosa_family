#!/bin/bash
# Author: Tao Xiong
# Date: 2025-12-16
# Description: 
# ==== 主体代码开始 ====

family=Meliaceae
outgroup_id=101569 
outgroup_name=Simarouba_amara

mkdir ${family} && cd ${family}

# Setup clade 可以使用parallel批量进行
python3 /home/miaosun/software/PyPHLAWD/src/setup_clade_ap.py -b /data/data1/Plant_DB2024.08/pln.202408.pyphlawd.db -s /data/data1/Plant_DB2024.08/plngzseqs/ -t ${family} -o . -l ${family}.log 
# parallel -j 8 '
#   mkdir -p {} &&
#   (
#     path=$(ls -d ./{}_*) &&
#     python3 /home/miaosun/software/PyPHLAWD/src/setup_clade_ap.py \
#       -b /data/data1/Plant_DB2024.08/pln.202408.pyphlawd.db \
#       -s /data/data1/Plant_DB2024.08/plngzseqs/ \
#       -t {} \
#       -o . \
#       -l {}.log
#   )
# ' :::: family_namelist.txt
# Finding good clusters
path=$(ls -d ./${family}_*)

yes y|python3 /home/miaosun/software/PyPHLAWD/src/find_good_clusters_for_concat.py -b /data/data1/Plant_DB2024.08/pln.202408.pyphlawd.db -d ./${path}| tee pyphlawd.log

#批量
# parallel -j 8 '
#   cd {} &&
#   (
#     path=$(ls -d ./{}_*) &&
#     yes y|python3 /home/miaosun/software/PyPHLAWD/src/find_good_clusters_for_concat.py -b /data/data1/Plant_DB2024.08/pln.202408.pyphlawd.db -d ./${path}| tee pyphlawd.log &&
#     cd /home/xiongtao/data/project/Loeess_Poisnoous/pyphlawd 
#   )
# ' :::: family_namelist.txt


# Finding min overlap
cmd=$(grep "^python3 .*get_min_overlap_multiple_seqs.py" pyphlawd.log)
OUTALN=$(ls ./${family}_*/|grep "outaln$")

# pxrmt pxrms
echo ./${family}_*/"$OUTALN" | eval "$cmd"|tee pxrmt.txt
cmd2=$(sed -z 's/.*pxrmt/pxrmt/' pxrmt.txt)
eval "$cmd2"

# Add outgroup
python3 /home/miaosun/software/PyPHLAWD/src/add_outgroup_to_matrix.py -b /data/data1/Plant_DB2024.08/pln.202408.pyphlawd.db -m ./${path}/${family}_*_outaln.filt -p ./${path}/${family}_*_outpart -t ${outgroup_id} -o ./${path}/${family}_OG -s /data/data1/Plant_DB2024.08/plngzseqs -mn 2 -mx 5

# rename -------------------------------------------------------------------------------------------
# Get the tips id from "${family}_26468.table"
cd ${path}
cut -f2,5 ${family}_[0-9]*.table|sort|uniq|cut -f1 >id.txt

#Get the tips species name
cut -f2,5 ${family}_[0-9]*.table|sort|uniq|cut -f2|sed 's/ /_/g'>sp_name.txt

# update the id and species name file by adding in the outgroup
echo $outgroup_id >>id.txt
echo $outgroup_name >>sp_name.txt

#Rename the seqence file
pxrls -s ${family}_OG.outaln -c id.txt -n sp_name.txt>${family}_OG_rn.outaln

#Rename the constrain tree
pxrlt -t ${family}_[0-9]*_outaln.constraint.tre.filt -c id.txt -n sp_name.txt>${family}_outaln.constraint_rn.tre.filt

mkdir final_tree && mv ${family}_outaln.constraint_rn.tre.filt final_tree/ && mv ${family}_OG_rn.outaln final_tree

cd /home/xiongtao/data/project/Loeess_Poisnoous/pyphlawd

# make tree -------------------------------------------------------------------------------------------
#将所有的的序列文件和限制树文件移动到一个文件夹中
find . -type f -path "*/final_tree/*" -exec mv -t all_tree {} +

#比对建树
mkdir mafft trimal genetrees
for name in ./*.outaln;do
  echo $(basename $name .outaln).mafft
  mafft --auto $name > mafft/$(basename $name .fasta).mafft
done

for name in ls mafft/*.mafft; do
  echo ${name}.tri.fasta
  output_file="trimal/$(basename "$name" .mafft).tri.fasta"
  trimal -in "$name" -out "$output_file" -automated1
done

for name in trimal/*.fasta;do
  tt=$(basename $name _OG_rn.tri.fasta)
  iqtree -s $name -g ${name}_outaln.constraint_rn.tre.filt -m MFP -B 1000 --bnni -T 10 -pre ./genetrees/${tt}
done