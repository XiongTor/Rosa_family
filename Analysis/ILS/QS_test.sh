#!/bin/bash
# Author: Tao Xiong
# Date: 2025-09-26
# Description: For QS test
# ==== 主体代码开始 ====
#用到的脚本位置：/data/xiongtao/scripts/quartetsampling/pysrc
#无需安装，在已安装python3和RAxML的基础上，脚本可直接运行
#quartet_sampling.py脚本用于进行Quartet Sampling（QS）计算。
#运行命令,记得指定路径
seq=$1
tree=$2
path=$(pwd)

mkdir -p ./{QS_result,QS_visual} && cd ./QS_result

python3 /data/xiongtao/scripts/quartetsampling/pysrc/quartet_sampling.py --tree ../$tree --align ../$seq -partitions ../$partitions --reps 1000 --threads 4 --lnlike 2 --verbout
   # –tree:    必需，指定拓扑树
   # -align:   指定的aligned文件
   # –reps:    每个枝计算的重复次数，默认100，推荐1000
   # –threads: 每个枝的重复的平行运算的线程，默认1
   # –lnlike:  log-likelihood阈值，默认2。用来指定最好似然树超过第二好似然树的最小差异（对每个枝来说）。推荐指定，如果不指定，则调用简单模式，用RAxML从比对序列推断树，不评估似然值，所有枝上的QI分数为0
   # –verbout: 指定后提供每个topology和QC的频率的输出，生成RESULT.verbout文件

Rscript /data/xiongtao/scripts/Topology_test/QS_visual.R $path