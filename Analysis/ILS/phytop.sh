#!/bin/bash
# Author: Tao Xiong
# Date: 2025-10-17
# Description: 繁
# ==== 主体代码开始 ====

# running ASTRAL use `-u 2`
astral-pro3 -i rosa_orthofinder_genetrees.tre -u 2 -o rosa_orthofinder_treeshrink_sp_u2.tre

# run phytop
phytop rosa_ags353_treeshrink_sp_u2.tre -figsize 1 -fontsize 4