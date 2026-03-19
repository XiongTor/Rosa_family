#!/bin/bash
# Author: Tao Xiong
# Date: 2026-03-19
# Description: 按照亚科，族水平倒塌树
# ==== 主体代码开始 ====
paste old_spname.txt Subfamily_name.txt >new_subf_label.map

nw_rename rosa_orthofinder_connect.txt.rt.subf.tre new_tribe_label.map >test.tre


for tree in sptree/*.tre;do
  tt=$(basename $tree .tre)
  nw_rename $tree new_tribe_label.map >test.tre
  nw_condense test.tre >tribe_tree/${tt}.tribe.tre
done