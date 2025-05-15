#!/bin/bash
# Author: Tao Xiong
# Date: 2025-05-15
# Description: This script is used to reroot the gene trees
# ==== 主体代码开始 ====

for name in *.treefile;do
  nn=$(basename $name .treefile)
  tt=$(nw_labels -I "$name" | grep -E "Zelkova_schneideriana|Morus_alba|Elaeagnus_pungens_hangzhou" | sed ':a;N;$!ba;s/\n/,/g')
  pxrr -t $name -g $tt>../treefile_rt/${nn}.rt.treefile
done
