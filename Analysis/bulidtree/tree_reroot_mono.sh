#!/bin/bash
# Author: Tao Xiong
# Date: 2026-03-19
# Description: 用于剪完树后的置根和单系性检测
# ==== 主体代码开始 ====

pxrr -t rosa_sp_partition.txt.treefile -g Elaeagnus_angustifolia,Zelkova_schneideriana,Morus_indica >rosa_sp_partition.txt.rt.tre

pxrlt -t rosa_sp_partition.txt.rt.tre -c old_name.txt -n new_name.txt > rosa_sp_partition.txt.rn.rt.tre

