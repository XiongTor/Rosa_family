#!/bin/bash
# Author: Tao Xiong
# Date: 2025-08-27
# Description: Genetree rerooting script
# ==== 主体代码开始 ====
for tree in ../auto/*.treefile; do
    nn=$(nw_labels -I $tree|grep -c -E "Elaeagnus_angustifolia_Armenia|SRR28089552|SRR27599619")
    if [[ $nn -gt 0 ]]; then
        tt=$(nw_labels -I $tree|grep -E "Elaeagnus_angustifolia_Armenia|SRR28089552|SRR27599619"|paste -sd,)
        mm=$(basename $tree .treefile)
        pxrr -t $tree -g $tt >${mm}.rt.tre
    else
        Rscript /home/xiongtao/data/scripts/mad_reroot.R $tree
    fi
done