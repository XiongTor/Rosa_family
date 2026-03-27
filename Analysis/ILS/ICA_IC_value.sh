#!/bin/bash
# Author: XiongTao
# Date: 2026-03-09
# Description: 
# ==== 主体代码开始 ====
/usr/bin/raxmlHPC-PTHREADS-AVX -T 20 -f i -t rosa_orthofinder_MO_treeshrink_sp_rt.tre -z rosa_orthofinder_MO_treeshrink_genetrees.tre -m GTRCAT -n T4 -C 