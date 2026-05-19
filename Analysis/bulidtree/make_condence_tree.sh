#!/bin/bash
# Author: Tao Xiong
# Date: 2026-03-19
# Description: 按照亚科，族水平倒塌树
# ==== 主体代码开始 ====
paste old_spname.txt Subfamily_name.txt >new_subf_label.map

nw_rename Rosaceae_MY_CP.tre new_genus_label.map >genus.tre


for tree in sptree/*.tre;do
  tt=$(basename $tree .tre)
  nw_rename $tree new_tribe_label.map >test.tre
  nw_condense test.tre >tribe_tree/${tt}.tribe.tre
done