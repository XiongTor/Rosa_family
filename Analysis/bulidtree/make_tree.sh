#各种方法建树,物种树

#modeltest
modeltest-ng -i rosa_sp_supermatrix.fasta -d nt -o modeltest

#树的拓扑结构对比 
cophylo(tre1=,tre2=,rotate=TRUE)

#raxml-ng
raxml-ng --all --msa rosa_sp_supermatrix.fasta --model GTR+I+G4 --outgroup Barbeya_oleoides --bs-trees 100 --threads 20
#先跑
raxml-ng --parse --msa ALIGNMENT --model GTR+G --prefix $Title --threads 1
#得到最佳模型mmmm，然后跑下面的
raxml-ng --all --msa ALIGNMENT.raxml.rba --model mmmm --tip-inner --prefix wwww --threads xxxx -seed ${RANDOM} --tree rand{25},pars{25} --bs-metric fbp,tbe


#iqtree
#GTR+F+R5,iqtree自选BIC
iqtree -s rosa_sp_supermatrix.fasta -m MFP -o Barbeya_oleoides -B 1000 --bnni -T AUTO
#GTR+I+G4，modeltest,AICc
iqtree -s rosa_sp_supermatrix.fasta -m GTR+I+G4 -o Barbeya_oleoides -B 1000 --bnni -T AUTO

#加分区文件的版本
iqtree -s Ela_sp_supermatrix.fasta --seed 12345 -o Barbeya_oleoides_100 -B 1000 -T 3 -p Ela_sp_partition.txt -m MFP -alrt 1000

#raxml
raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 10 -s rosa_sp_supermatrix.fasta -n rosa.raxml -T 15

#-f a 此参数用于选择 RAxML 运算的算法。可以设定的值非常之多。 a 表示执行快速 Bootstrap 分析并搜索最佳得分的 ML 树。

#-x 12345 指定一个 int 数作为随机种子，以启用快速 Bootstrap 算法。

#-p 12345 指定一个随机数作为 parsimony inferences 的种子。

#-# 10 指定 bootstrap 的次数。

#-m GTRGAMMA 指定核苷酸或氨基酸替代模型。 "PROT" 表示氨基酸替代模型，“GTR”表示碱基替代模型； GAMMA 表示使用 GAMMA 模型； X 表示使用最大似然法估计碱基频率。

#-s 20k.phy 指定输入文件。phy 格式的多序列比对结果。软件包中包含一个程序来将 fasta 格式转换为 phy 格式。也可以通过Tassel或者Mega转换格式：vcf-phylip

#-n chr001.raxml 输出文件的后缀为 .chr001.raxml 。

#-T 30 指定多线程运行的 CPUs 。


#ASTRAL
cat ./Tree_reroot/rosatree*> rosa_sp_genes.tre
#Astral建立物种树
java -jar /home/xiongtao/data/software/ASTRAL-master/Astral/astral.5.7.8.jar -i rosa_sp_genes.tre --outgroup Barbeya_oleoides -o rosa_Astral_species.tre

#重新估算枝长
raxmlHPC -f e -t rosa_Astral_species_sortadata.tre -m GTRGAMMA -s ../rosa_sp_supermatrix_sortadata.fasta -n rosa_astral_species_br.tre

#给物种树置根
pxrr -t rosa_sp_supermatrix.fasta.raxml.support -g Barbeya_oleoides >rosa_sp_supermatrix.fasta.raxml.rt.support



#Bayes
BEGIN MRBAYES;
outgroup Bolusia_amboensisMlR112 #外群名1
outgroup Lotononis_hirsuta #外群名2
lset nst=6 rates=invgamma ngammacat=4; #进化模型参数
Prset statefreqpr=dirichlet(1,1,1,1); #进化模型参数
mcmc ngen=16000000 printfreq=1000 nruns=2 diagnfreq=5000 samplefreq=10000 nchains=4 temp=0.1 burninfrac=0.25 checkpoint=yes savebrlens=yes;
sump burnin=400;
sumt burnin=400 contype=allcompat showtreeprobs=yes;
END;

