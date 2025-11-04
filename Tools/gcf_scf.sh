#!/bin/bash
# Author: Tao Xiong
# Date: 2025-10-23
# Description: For count gcf and scf value
# ==== 主体代码开始 ====
iqtree3 -t rosa_ags353_treeshrink_sp_rt.tre --gcf rosa_ags353_treeshrink_genetrees.tre -p ../02-trimal --scf 100 --prefix concord -T 15
