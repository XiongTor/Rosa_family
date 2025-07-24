#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-05-15
# Description: This script is used to run the PhyloNetworks (SnaQ) analysis

# ==== Main Script Start ====

# base julia to process the data

# install julia

# 1.use linux pkg
wget https://julialang-s3.julialang.org/bin/linux/x64/1.7/julia-1.7.2-linux-x86_64.tar.gz #下载linux的64位预编译的julia
tar -xzf julia-1.7.2-linux-x86_64.tar.gz #解压缩
julia-1.7.2/bin/julia -h #查看帮助文章

# 2. use julia github code
git clone https://github.com/JuliaLang/julia.git #克隆github上的最新版源代码
cd julia
git checkout v1.7.2 #运行checkout来获取julia的最新稳定版本1.7.2
make #编译
./julia -h #查看帮助文档


#################### use julia beginning #################################
# install PhyloNetworks pkg
using Pkg-Pkg.add("PhyloPlots")

# prepare the input data
# 1.CF(consistency factor) table, means Gene Tree Frequency Table
# 2. start tree ,you can use the species tree

####################### use gene trees to get the CF table#####################
#导入包
using PhyloNetworks 
using PhyloPlots 
using CSV 

#如果PhylpPlots包报错是R库的问题，请尝试配置环境变量，使得它能够找到正确的R库，注意在bash中运行
# export LD_LIBRARY_PATH="/home/xiongtao/miniconda3/envs/r440/lib:/home/xiongtao/miniconda3/envs/r440/lib/R/lib:$LD_LIBRARY_PATH"


#读取基因树文件
raxmltrees=joinpath("rosa_tribe_rt_genetrees.tre") 

#解析基因树
genetrees = readMultiTopology(raxmltrees) 

#读取基因树，计算四分类群的CFs
q,t = countquartetsintrees(genetrees)

#读取计算得到的CF值到df：基因频率
df = writeTableCF(q,t) 

 #保存df内容为tableCF.csv文件
CSV.write("tableCF.csv", df)

######################################## make the network tree ###############################################
#导入起始树文件
astralfile= joinpath("rosa_tribe_Astral_species.rt.tre") 

# 读取树文件
astraltree = readMultiTopology(astralfile)[1]

# 读取tableCF.csv文件，生成"DataCF" 对象
raxmlCF = readTableCF("tableCF.csv") 

# 用astraltree,raxmlCF，指定hmax=0来运行。
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234) 

#直接运行接下来的若干net即可，可以尝试使用脚本并行或者迭代运行
