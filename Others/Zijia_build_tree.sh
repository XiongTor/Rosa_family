## 将下载的序列整合为基因-物种形式
#示例脚本：sun_organized_seq.sh
#!/bin/bash
#author: Miao Sun
# this script will parse equences from one of 353 genes based on each species
# then distribute alignment of each gene under each genus folder 
# then serve for the phylogeny for each gene for each genus
# using mafft to aligment; if the some sequences is reverse direction, then 
# automatically reverse-complemented it

#seqlist=$1
#genelist=$2

if [ ! -e "reformat_aln" ]; then #若当前目录没有”reformat_aln" 文件夹，则创建该目录用于存储结果
        mkdir reformat_aln
fi

while read -r Line; do  #读取species_list.txt每个文件路径

        seqfile=$(basename $Line .fasta)  #根据文件路径提取文件名，不包含路径和fasta后缀

        echo -e "\n\n$seqfile\n\n"
        
        while read -r gene; do
        
                # check if a gene sequence file exists
                if [ ! -f ./reformat_aln/${gene}.fasta ]; then
                        touch ./reformat_aln/${gene}.fasta #为每个基因创建空输出文件
                fi
                
                II=$(grep $gene $Line|wc -l) #统计当前文件中含有该基因ID的条目数量
                
                if [[ "$II" -gt 0 ]]; then
                         grep -EA 1 ">$gene" $Line|sed -E 's/>.* Species:([^ ]+).* Sequence_ID:([^ ]+).*/>\1_\2/' >> "./reformat_aln/${gene}.fasta" #提取基因ID所在行及下一行（序列行），只保留序列头中的species_ID，添加">$gene" 以索引到每个基因
                fi
        done <Ags353_gene_list.txt  #每行一个基因名，对应原始文件序列头的基因名如：4471
#done<species_list.txt
done<species_list.txt #该文件为每行一个FSATA文件路径：./paftol_download/INSDC.DRR238827.Weddellina_squamulosa.a353.fasta

##==============================================================================================================================
# mafft，将每个353基因序列进行比对，保存为gene.mafft 文件
#示例脚本：mafft_auto.sh
mkdir mafft
for gene in ./reformat_aln/*.fasta
do
# Perform multiple sequence alignment and save output to file
    mafft --auto --thread 15 "$gene" > "${gene%.*}.mafft"
    mv "${gene%.*}.mafft" ./mafft
done

# prank 精细比对
#1.先conda 安装PRANK
#2.prank_wrapper.py
mkdir prank
conda activate PRANK
python /data/miaosun/script/PRANK/prank_wrapper.py ./mafft ./prank mafft dna

# 得到的结果prank 共349个基因，mafft 有353个，看看哪些基因被去掉了
ls ./mafft/*.fasta >mafft.txt
ls ./reformat/* >prank.txt
comm -3 \
  <(cut -d '.' -f1 prank_list | sort -u) \
  <(cut -d '.' -f1 trimal_auto_list | sort -u) \
| sed '/^\s*$/d; s/\t//g'
# 结果得到共有4个基因：5034、6164、6216、6221这几个基因没有在prank结果中

#统计每个基因文件有多少个物种序列
#示例脚本
for f in reformat_aln/*.fasta; do
    bn=$(basename $f .fasta)
    cnt=$(grep '^>' $f | cut -d'_' -f1-2 | sort -u | wc -l)
    echo "$bn $cnt"
done | sort -k2nr > species_count_sorted.txt

# 发现有一个基因6514只有一个物种：Weddellina_squamulosa
# 反向查找所有物种是不是只有这一个物种有这个基因
grep -h '>6514' paftol_download/*.fasta
# >6514 Gene_Name:PPOX2 Species:Weddellina_squamulosa Repository:INSDC Sequence_ID:DRR238827 结果就这一个物种

# 去掉这物种后，使用352个基因用于后续构建进化树，首先进行mafft、prank 比对、trimal修剪
#mafft 比对352条，prank比对352条，trimal修剪
# 使用了五种参数尝试修剪：
# trimal_auto1
mkdir ./trimal_auto1
for gene in ./prank/*.mafft.aln; do
    trimal -in "$gene" -out "${gene}.tri.fasta" -automated1
    mv "${gene}.tri.fasta" ./trimal_auto1
done

# trimal_gt0.8_cons60_st0.001
mkdir ./trimal_gt0.8_cons60_st0.001
for gene in ./prank/*.mafft.aln; do
    trimal -in "$gene" -out "${gene}.tri.fasta" -gt 0.8 -cons 60 -st 0.001
    mv "${gene}.tri.fasta" ./trimal_gt0.8_cons60_st0.001
done

# trimal_gt0.2cons85_st0.001
mkdir ./trimal_gt0.2cons85_st0.001
for gene in ./prank/*.mafft.aln; do
    trimal -in "$gene" -out "${gene}.tri.fasta" -gt 0.2 -cons 85 -st 0.001
    mv "${gene}.tri.fasta" ./trimal_gt0.2cons85_st0.001
done

# trimal_gappyout
mkdir ./trimal_gappyout
for gene in ./prank/*.mafft.aln; do
    trimal -in "$gene" -out "${gene}.tri.fasta" -gappyout
    mv "${gene}.tri.fasta" ./trimal_gappyout
done

# trimal_strictplus
mkdir ./trimal_strictplus
for gene in ./prank/*.mafft.aln; do
    trimal -in "$gene" -out "${gene}.tri.fasta" -strictplus
    mv "${gene}.tri.fasta" ./trimal_strictplus
done

#先用以上修剪结果建树，以及再用prank 结果建树
#=========================================================================================================================
#raxmlHPC建树，但是raxmlHPC模型较少，没有raxml-ng模型多
#!/bin/bash

# 创建一个目录用来存放生成的系统发育树
mkdir -p gene_tree_trimalauto1
cd gene_tree_trimalauto1

# 遍历 trimal 目录下的每个 .fasta 文件
for file in ../trimal_auto1/*.fasta; do
echo "Processing $(basename "$file")"

# 使用 RAxML 构建系统发育树并命名输出文件
raxmlHPC -f a -s "$file" -k -x $RANDOM -m GTRGAMMA -p $RANDOM -n "$(basename "$file" .fasta).ML.tre" -w "$(realpath .)" -N 100
done
#=============================重要注意===========================================================================
#建树时最好不要直接指定GTRGAMMA，应该用iqtree或者raxml-ng自己指定最佳模型
#其次使用raxml-ng，不使用raxml-ng-mpi
#modeltest
mkdir -p modeltest_results
cd modeltest_results
for file in ../trimal_auto1/*.fasta; do
    prefix=$(basename "$file" .fasta)
    echo "==== ModelTest for: $prefix ===="
    
    # 用 ModelTest-NG 选择最优模型
    modeltest-ng -i "$file" -d nt -o "$prefix" -p 1 
done
# 提取每个模型测试结果
echo "gene best_model">./model_list.txt

for name in *.out;do
    echo $name
     tt=$(basename $name .mafft.aln.tri.out)
     nn=$(grep -A 2 "Best model according to AICc" $name|tail -n 1|cut -f2 -d':'|sed 's/\s//g')
    echo -e "$tt $nn">>./model_list.txt
done

#raxml-ng模型测试：

mkdir -p raxml-ng_model_test
cd raxml-ng_model_test

for file in ../trimal_auto1/*.fasta; do
    echo "Processing $(basename "$file")"

    raxml-ng --parse --msa "$file" --model GTR+G --prefix ./$(basename "$file".fasta) --threads 1 

done
# 提取模型结果
echo "gene model">./model_list.txt

for name in *.log;do
    echo $name
     tt=$(basename $name )
     nn=$(grep "Model:" *.log  |cut -f2 -d' ')
    echo -e "$tt $nn">>./model_list.txt
done
#raxml-ng 模型检测
raxml-ng --parse --msa ./trimal/*.fasta --model GTR+G --prefix $Title --threads 1

#raxml-ng 建树
mkdir -p raxml-ng_trimal_auto1
cd raxml-ng_trimal_auto1

    # 遍历 trimal 目录下的每个 .fasta 文件
    for file in ../trimal_auto1/*.fasta; do
    echo "Processing $(basename "$file")"

    # 使用 RAxML 构建系统发育树并命名输出文件
    raxml-ng --msa "$file" --model GTR+FO+G4m --prefix ./$(basename "$file".fasta) --threads 16 --tip-inner on --seed ${RANDOM} --tree rand{25} --bs-metric fbp
    
done
#直接使用iqtree也是方便快捷准确的选择====================================================================================================================
#脚本示例：iqtree_trimal_auto1.sh
# Iqtree
mkdir -p Iqtree_gene_tree_trimalauto1/
cd Iqtree_gene_tree_trimalauto1/

# 遍历 trimal_auto1 目录下的每个 .fasta 文件
for file in ../trimal_auto1/*.fasta; do
    prefix=$(basename "$file" .fasta)  
    echo "Processing $prefix"
    
    #使用iqtree选择最优模型并构建系统发育树
    iqtree2 -s "$file" -m MFP --alrt 1000 -B 1000 -T 10

done
 
#=========================================================================================================================
## 收集并整理基因树分析结果，将最佳进化树文件(*.treefile) 及其对应的 FASTA 比对文件归类到特定目录，便于后续分析（如构建物种树、数据集整合等）
#脚本示例：clear_up_genetree_genealn.sh
#!/bin/bash

# 创建目标目录
mkdir -p Iqtree_trimalauto1_best_tree Iqtree_trimalauto1_best_aln

# 遍历所有基因树文件（.treefile）
for treefile in Iqtree_trimal_auto1/*.treefile; do
    # 提取基因名称（例如：7628.mafft.aln... → 7628）
    gene_name=$(basename "$treefile" | grep -o '^[0-9]\+')
    
    # 1. 复制基因树文件到 gene_best_tree 并简化名称
    cp -v "$treefile" "Iqtree_trimalauto1_best_tree/${gene_name}.treefile"
    
    # 2. 定义对应的 FASTA 文件路径
    fasta_src="trimal_auto1/${gene_name}.mafft.aln.tri.fasta"
    fasta_dest="Iqtree_trimalauto1_best_aln/${gene_name}.fasta"
    
    # 3. 检查并复制 FASTA 文件
    if [[ -f "$fasta_src" ]]; then
        cp -v "$fasta_src" "$fasta_dest"
    else
        echo "警告: 未找到文件 $fasta_src"
    fi
done

echo "整理完成！结果文件保存在 gene_best_tree 和 gene_aln_best 目录中。"

#=========================================================================================================================
#置根每个基因树
#整体思路：外类群与内类群均不为单系，无法置根，因此下面置根策略都是选择一个外类群置根，只不过制定了优先级别：
#如果在这棵树中四个外类群都有，选择Gunnera_manicata_SRR7451106置根，如果没有这个物种，其次选择是Dillenia_indica_EHNF、Penthorum_chinense_ERR4180013、Sarcophyte_sanguinea_ERR7621663。
#最后四个物种都没有的基因树选择mad法置根
# 脚本在Reroot_Gene_Trees.sh
#!/bin/bash

# 1. 创建输出目录
mkdir -p Tree_reroot4

# 2. 遍历所有基因树文件（.treefile）
for tree in Iqtree_trimalauto1_best_tree/*.treefile; do 
    echo "正在处理文件: $tree"

    # 3. 提取基因树名称
    gene_name=$(basename "$tree" .treefile) 
    treename="${gene_name}"  # 输出文件前缀

    # 4. 检查树中包含的外类群
    contains_GM=$(pxlstr -t "$tree" -i | grep -q "Gunnera_manicata_SRR7451106" && echo 1 || echo 0)
    contains_DI=$(pxlstr -t "$tree" -i | grep -q "Dillenia_indica_EHNF" && echo 1 || echo 0)
    contains_PC=$(pxlstr -t "$tree" -i | grep -q "Penthorum_chinense_ERR4180013" && echo 1 || echo 0)
    contains_SS=$(pxlstr -t "$tree" -i | grep -q "Sarcophyte_sanguinea_ERR7621663" && echo 1 || echo 0)
    
    # 5. 选择优先级最高的外群
    if [[ "$contains_GM" -eq 1 ]]; then
        OG="Gunnera_manicata_SRR7451106"
    elif [[ "$contains_DI" -eq 1 ]]; then
        OG="Dillenia_indica_EHNF"
    elif [[ "$contains_PC" -eq 1 ]]; then
        OG="Penthorum_chinense_ERR4180013"
    elif [[ "$contains_SS" -eq 1 ]]; then
        OG="Sarcophyte_sanguinea_ERR7621663"
    else
        OG=""
    fi

    # 6. 置根操作
    if [[ -n "$OG" ]]; then
        echo "使用 $OG 进行置根..."
        pxrr -t "$tree" -g "$OG" -o "./Tree_reroot4/${treename}.rt.tre"
    else
        echo "未检测到外群，使用 MAD 方法重根..."
        Rscript ./mad_reroot.R "$tree"
    fi

done

#以上脚本已验证，所有基因树都成功置根
#统计置根情况：if_genetree_rooted.sh
for file in `ls *.xxx`; 
    do nn=$( -t $tree -n);
    rr=$(pxlstr -t $tree | grep "rooted" | cut -f 2 -d " ");
    echo -e "$tree\t$rr\t$nn" ; 
done
# 以上构建好所有的基因树了，可以通过nw_display 看进化树的结构
nw_display .treefile
#======================================================================================================================================
#做完基因树后有许多树有长枝，treeshrink 用来去除一些异常值
#安装
git clone https://github.com/uym2/TreeShrink.git

python setup.py install

#尝试是否安装成功，查看参数
run_treeshrink.py -h

# **1.序列比对文件和系统发育树的准备**
# treeshrink_input_prepare.sh
#!/bin/bash

# 遍历所有树文件(.rt.tre)并处理
for tre_file in reroot_genetree_trimalauto1_iqtree/*.rt.tre; do
    # 提取基因 ID（去掉 ".rt.tre" 后缀）
    gene_id=$(basename "$tre_file" .rt.tre)
    fasta_file="best_aln_iqtree_trimalauto1/${gene_id}.fasta"
    
    # 检查对应的 fasta 文件是否存在
    if [ -f "$fasta_file" ]; then
        # 创建基因文件夹
        mkdir -p "$gene_id"
        
        # 复制文件到目标目录（保留原文件）
        cp "$tre_file" "$gene_id/input.tree"      # 树文件复制为 input.tree
        cp "$fasta_file" "$gene_id/input.fasta"   # 比对文件复制为 input.fasta
        echo "已处理基因: $gene_id"
    else
        echo "[错误] 文件缺失: $fasta_file 未找到"
    fi
done
#文件夹结构：
#输入
# 当前目录/
├── reroot_genetree_trimalauto1_iqtree/
│   ├── 7628.rt.tre    # 原始未修改的树文件
│   ├── 6864.rt.tre
│   └── ...
└── best_aln_iqtree_trimalauto1/
    ├── 7628.fasta     # 原始未修改的比对文件
    ├── 6864.fasta
    └── ...
# 输出
├── 7628/                                # 新增基因子目录
    ├── input.tree    # 来源于 reroot_genetree_trimalauto1_iqtree/7628.rt.tre 的副本
    └── input.fasta   # 来源于 best_aln_iqtree_trimalauto1/7628.fasta 的副本
├── 6864/
    ├── input.tree
    └── input.fasta
...

#Treeshrink检测异常值仅仅是基于input tree，序列不影响异常值的检测，github网址：https://github.com/uym2/TreeShrink
#需要提供一个文件夹，如gene1, gene2,...,gene10，文件夹内包括树和矩阵：input.tree和input.fasta。文件夹内部的树和矩阵名称相同
#-t gene tree
#-a alignment files
# -m 三种参数：
# per-gene'：这种模式仅在输入的树之间是系统发育独立时才有用
# 'all-genes'
# 'per-species'：优先级最高的是 'per-species' 模式，除非数据集中基因树的数量太少（即少于 20 棵树）或存在稀有物种（即某个物种在少于 20 棵树中出现）
# -b 若某一物种对直径的影响小于给定值 x%（x 为数值）时，在基于单个物种的测试中则不将其移除。默认值：5 。github建议探索20
#-i 在放置所有基因文件夹的目录下即可，我的是/home/zijia/data/mk_tree/download_fasta/treeshrink/reroot_trimalauto1_iqtree

nohup run_treeshrink.py -i ./ -t input.tree -a input.fasta -m per-species -b 20 > ./input.tree.treeshrinklog.txt &
#结果在每个gene/文件夹下生成了output.fasta（建树的基因序列，去除了物种）、output.tree（去除物种后的树）、output.txt（检测到异常值去除的物种）
#output.tree依然是有根的
#整理结果文件
#!/bin/bash
mkdir -p treeshrink_aln treeshrink_tree

for dir in */                           # 遍历所有目录
do 
    gene=${dir%/}                       # 去除目录末尾的"/"（如将 "gene/ → gene）
    cp "$dir/output.fasta" "treeshrink_aln/$gene.fasta"
    cp "$dir/output.tree" "treeshrink_tree/$gene.tree"
done

#==========================================================================================================================================
# zuntini参数
python $APPS/TreeShrink/run_treeshrink.py \
-t ${PREFIX}_$iteration.trees \
--centroid \ # 在预处理过程中启用质心重根化。对于大型树种，强烈推荐使用此功能。默认设置：关闭
-m 'all-genes'
#=================================================================================================================================================
# cat合并树，倒塌10
cat ./reroot_genetree_trimalauto1_iqtree/*.tre> ./cat_genetree/cat_trimalauto1_gene.tre
#低于10的倒塌
nw_ed cat_trimalauto1_gene.tre 'i & b<=10' o > trimalauto1_BS10.treefile
# 为什么昨晚treeshrink 倒塌10 所有的树都变成梳子，因为这棵树没有支持率
# astral 合并物种树，（这里outgroup怎么选择）
java -jar /home/zijia/data/software/ASTRAL-master/Astral/astral.5.7.8.jar -i trimalauto1_BS10.treefile --outgroup Gunnera_manicata_SRR7451106 -o trimalauto1_BS10_Astral.tre
#置根
nw_reroot trimalauto1_BS10_Astral.tre Gunnera_manicata_SRR7451106>trimalauto1_BS10_Astral.rt.tre

#合并构建基因树用的所有基因为supermatrix，（这里supermatrix是treeshrink后生成的新fasta（对，如果做了treeshrink的话），cat之后需要重新生成新fasta吗（不需要，因为只是把枝收起来了））
mkdir -p supermatrix
pxcat -s Iqtree_trimalauto1_best_aln/*.fasta -p Iqtree_trimalauto1_best_aln_partition.txt -o ./supermatrix/Iqtree_trimalauto1_best_aln_supermatrix.fasta
#生成partition.txt、supermatrix.fasta
#重新估算枝长
iqtree2 -s ./supermatrix/Iqtree_trimalauto1_best_aln_supermatrix.fasta -p ./supermatrix/Iqtree_trimalauto1_best_aln_partition.txt -m MFP -g ./astral_species_tree/trimalauto1_BS10_Astral.tre --prefix trimalauto1_BS10_Astral_branch_length >trimalauto1_BS10_Astral_branch_length.txt 2>&1

#置根
nw_reroot trimalauto1_BS10_Astral_branch_length.treefile Gunnera_manicata_SRR7451106 >trimalauto1_BS10_Astral_branch_length_re.tre

#最后结果为xxBS10_Astral_branch_length_re.tre