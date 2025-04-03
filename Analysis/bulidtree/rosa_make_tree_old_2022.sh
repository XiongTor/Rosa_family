ar#!/bin/bash
# Author:xiongtao
#data : 2023.05.16
# what function to use "ML" (dating the best ML tree) or "BS" (dating the bootstrap tree)
action=$1
#本流程用到的所有脚本，均在本地的  《 D:\R\代码学习脚本\笔记\蔷薇科建树脚本 》 中有备份。
#1.数据来源与处理：
#主要来源于paftol下载的数据与easy353自己筛选的数据：
#需要将paftol中一个物种的每一个353基因单独摘出来，放到一个新的文件中，为此写了一个脚本，放到本地的“蔷薇科建树脚本”文件夹中的“paftol基因摘出.sh”。
#服务器位置如下：/data/xiongtao/scripts/paftol_353gene.sh
#具体脚本内容如下：
#================================================================================
#将所有fasta文件的多行序列转化为单行序列,并去掉为空行的第一行,同时将原paftol文件和转化后的paftol文件移动到新的文件夹中。

#创造新的文件夹
    mkdir basedata && cd basedata && mv ../* ./
    mkdir paftol_new

    #去掉paftol文件中的换行符
    for name in *;do
       if [ -f "$name" ]; then
          awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' $name |sed 1d - > ${name}.fasta
          mv $name paftol_new
       fi
    done

    #摘出各个物种的基因名
    for file in ./*.fasta;do 
      gg=$(basename $file .fasta)
      grep ">" $file | awk -F'[ :]' '{print $1}'|sed 's/>//g'>${gg}_gene_list.txt

    #创建目录用于存储分离出的基因文件
      mkdir ${gg}
      echo ${gg}

    #逐行读取基因列表文件
        while read -r gene;do
    #使用grep命令查找基因名称，并将匹配的结果保存到单独的文件中
            grep -A 1 -w "${gene}" "$file"> "${gg}/$(echo ${gene}).fasta";
        done < "${gg}_gene_list.txt"
    done

    mv *.fasta paftol_new
    rm *.txt

    #之后需要确认easy353筛选出来的文件和paftol文件的区别，确保他们的一致性。
    #本次流程以paftol文件为基准，调整easy353文件，例如将easy353文件中的target文件夹提取出来命名为物种名，
    #将target文件中的target.fasta改为.fasta。总之本步后续看各人文件不同，自行调整为一种类型即可。

    #针对easy353筛选出来的文件，将target_genes中的文件放到一个新的，用物种名命名的文件夹中。
    for directory in *;do
        if [[ ! -d "./$directory/target_genes" ]]; then
            echo paftol
        else
            rm -rf $directory/filtered_reads
            echo $directory &&mkdir ${directory}_1
            mv ${directory}/target_genes/*.fasta ./${directory}_1
            rm -rf ${directory}
            mv ${directory}_1 ${directory}
        fi
    done

    #统一基因文件的名称，去掉easy353文件中的target
    for file in $(ls -d */);do 
      cd ${file}
        for name in *.fasta; do 
          if [ -f *.target.fasta ]; then
	        mv $name $(echo $name|sed 's/.target//g');
          else
            echo "nothing"
          fi
        done
      cd ../
    done

    #同时需要调整paftol和353筛选出来的fasta文件的标题行，保证只有物种名和序列号，如：>Purshia_tridentata

    #easy353
    for file in $(ls -d */);do
    cd $file
    for name in *.fasta;do 
    #easy353
    sed -i 's/_[0-9]*_recovered.*$//g;s/_[0-9]\.[0-9]//g' $name
    #paftol
    #sed -i -E 's/>([0-9]+) Gene_Name:[^ ]+ Species:([^ ]+) Repository:[^ ]+ Sequence_ID:[^ ]+/>\2/' $name
    #sed -i 's/^\(>[a-zA-Z]\+\)_\([a-zA-Z]\+\).*/\1_\2/' $name
    #mv $name g$name
    done
    done
    cd ..



#============================================================================================================



#============================================================================================================

#调整完毕后即可将全部353基因按基因名合并为新的fasta文件。脚本放到本地，命名为  conbining_easy353_gene.sh。位置在：/data/xiongtao/scripts/conbining_easy353_gene.sh.具体内容如下:
#===========================================================================================
#!/usr/bin/bash
# Author:xiongtao
# This script is used for extracting each specific gene from each species, and combining them as one gene alignment
    mkdir conbing_gene && cd conbing_gene
# workdirectory
#Wpath=`pwd`
#Wpath="/home/xueqin/data/Prunus_dada_analyses/Transcriptome"

#cd $Wpath 

#create a 353gene name list
    ls /data/xiongtao/easy353/Easy353/353_ref_rosaceae/353gene > Ags353_gene_list.txt

   # output folder
    mkdir -p gene_alns && cd ./gene_alns/ && rm -rf *.fasta

   # creat fasta file name for each gene 
   # and fill in them later 
    while read -r Line; do
       touch g$Line; 
    done<../Ags353_gene_list.txt

   # create a species list to work on 
    ls ../../basedata|grep "_" >../species_list.txt
   # delete none species names
    sed '/Ags353_gene_list.txt/d' ../species_list.txt > ../species_name_list.txt


    rm ../species_list.txt

   # loop through the whole species list to extract genes
    while read -r species; do
   #species=$1
    echo $species
   #species="Maddenia_wilsonii"
    cd ../$species
      for gg in `ls *.fasta`; do
         Gene=$(echo $gg|sed 's/.fasta//g')
         echo $Gene
         cat $gg >>../gene_alns/${Gene}.fasta
      done
    done <../species_name_list.txt

    cd ../..
   #如果第31行报错，请检查是不是species_list.txt中有空行或者名称错误
#==============================================================================================================

#安装mafft进行排序与比对
#教程： https://www.jianshu.com/p/2b882849ed61
#循环比对
  mkdir mafft trimal
  for name in `ls gene_outlier`;do
    echo $(basename $name .fasta).mafft
    mafft --auto gene_outlier/$name > mafft/$(basename $name .fasta).mafft
  done

#循环过滤
  for name in ls mafft/*.mafft; do
    echo ${name}.tri.fasta
    output_file="trimal/$(basename "$name" .mafft).tri.fasta"
    trimal -in "$name" -out "$output_file" -automated1
  done

#=============================================================================================================

#循环建树并置根
    mkdir gene_tree && cd gene_tree 
    for file in ../continue/*.fasta;do
      echo $(basename "$file" trimal/)
      raxmlHPC -f a -s $file -k -x $RANDOM -m GTRGAMMA -p $RANDOM -n $(basename "$file" trimal/).ML.tre  -w $(realpath .) -N 100
    done

    #批量进行基因树重命名
    mkdir gene_best_tree gene_aln_best
    cp -r gene_tree/RAxML_bestTree* gene_best_tree
    #将所有基因树的名称摘出并重命名，保持和需要用来进行数量比对的fasta基因矩阵文件的名称一致。例如：6000.fasta
    ls gene_best_tree > besttree.txt
    sed "s/RAxML_bestTree.//g;s/.ML.tre//g" besttree.txt >besttree_new.txt
    #到存放fasta基因矩阵的文件目录下，将能成功比对上名字的文件复制出来到一个新的文件夹中
    parallel cp trimal/{} gene_aln_best/{} <besttree_new.txt
    #到此，所有能比对上名字的fasta基因矩阵都在rosa_gene_aln_best目录下，后续可以在这个目录下合成矩阵用于物种树建立

    # make a folder for rooted trees
    mkdir -p Tree_reroot

    # loop through each tree if the outgroup "Barbeya_oleoides"
    # prsent the use `pxrr` or `nw_reroot`
    for tree in `ls gene_best_tree/RAxML_bestTree*`; do 
      echo $tree
      # define a base name for rooted tree
      treename=$(basename $tree .tri.fasta.ML.tre|sed 's/RAxML_bestTree./rosatree./g')

      #echo $treename; 
      # If the outgroup "Barbeya_oleoides" prsent in the tree
      # then NN=1,else =0
      NN=$(pxlstr -t $tree -i|grep -c "Barbeya_oleoides")

      # creat in if condition

    # if outgroup present
      if [[ "$NN" -eq 1 ]]; then
        #echo -e "\n*******\nThis tree has outgroup using nw_reroot to root\n*******\n"
        #nw_reroot $tree Barbeya_oleoides >./Rosa_tree_reroot/${treename}.rt.tre
        # get outgroup from each tree
        OG=$(pxlstr -t $tree -i|grep "Barbeya_oleoides")

        # reroot using phx
        pxrr -t $tree -g $OG -o ./Tree_reroot/${treename}.rt.tre
      else # if the outgroup absent
        echo -e "\n*******\nNo outgroup found for $tree using mad.R to root\n*******\n"
        Rscript /home/xiongtao/data/scripts/mad_reroot.R $tree
      fi
  done

#接着上述提到的第一
#合成矩阵，注意，要求比对后的fasta文件的信息行/标题行中只有物种名，否则合成的矩阵中可能该物种的全部序列不会合到一起而是区分开了。
#然后可以开始合成矩阵了：
    mkdir species_tree
    pxcat -s gene_aln_best/*.fasta -p rosa_sp_partition.txt -o rosa_sp_supermatrix.fasta

    #gene_aln_best 包含有所有成功建立了基因树的比对清理后的fasta文件（即mafft和trimal之后的文件）
    #-p rosa_sp_partition.txt  该代码自动生成的分区文件
    #-o rosa_sp_supermatrix.fasta 生成的矩阵，包含所有物种的所有基因序列。（其中每个物种应该只出现一次）

    #合并置根之后的基因树
    cat ./Tree_reroot/rosatree*> rosa_sp_gene.tre
    #Astral建立物种树，Astral下载地址：https://github.com/smirarab/ASTRAL
    java -jar /home/xiongtao/data/software/ASTRAL-master/Astral/astral.5.7.8.jar -i rosa_sp_gene.tre --outgroup Barbeya_oleoides -o rosa_Astral_species.tre

    #重新估算枝长
    raxmlHPC -f e -t rosa_Astral_species_sortadata.tre -m GTRGAMMA -s ../rosa_sp_supermatrix_sortadata.fasta -n rosa_astral_species_br.tre

    #给物种树置根
    pxrr -t rosa_sp_supermatrix.fasta.raxml.support -g Barbeya_oleoides >rosa_sp_supermatrix.fasta.raxml.rt.support
#==================================================================================================================================================================

    #phyparts分析
    /home/xiongtao/data/tree/rosa_tree/rosa_best_gene_tree

    java -jar /data/xiongtao/software/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d sor_gene_tree -m Astral/RAxML_result.rosa_astral_species_br.rt.tre -o ./phyparts 

    bash /data/xiongtao/scripts/phyparts_pies.sh


#=============================================================================================================================================================


#iqtree构建系统发育树
iqtree -s NB-ARC.domain.afa  -nt 50
# -s 输入文件，-nt CPU 线程数

#查看当前目录占用的硬盘大小
du -sh

iqtree -s rosa_sp_supermatrix.fasta -m MFP -o Barbeya_oleoides -B 1000 --bnni -T 20


