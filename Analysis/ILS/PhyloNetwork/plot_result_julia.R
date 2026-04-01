#!/usr/bin/env julia
# Author: Tao Xiong
# Date: 2026-02-12
# Description: For displaying the results of PhyloNetworks analysis

# ==== Main Script Start ====
using PhyloNetworks
using PhyloPlots #注意更新到最新版julia,>1,9否则可能报错
using RCall

net = readTopology("net3_10runs.out")  # 将最佳网络树存为变量net
rootatnode!(net,"Morus_indica")  # 对树进行定根
writeTopology(net,di=true, "bestnet_h3.tre") # 将最佳网络树写到bestnet_h3.tre文件,并保存为dendroscopre 可读的形式

# 可以直接到dendroscopre中可视化，或者按照如下流程在julia中继续操作
imagefilename = "snaqplot_net3_root.svg"
R"svg"(imagefilename, width=8, height=6)
R"par"(mar=[0,0,0,0]) 
plot(net,
     showgamma=true,
     tipcex=0.9,
     xlim=[0,30])  # 横向空间扩大
R"dev.off()";

