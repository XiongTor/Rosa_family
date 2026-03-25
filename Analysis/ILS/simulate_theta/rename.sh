#!/bin/bash
# Author: Tao Xiong
# Date: 2026-03-23
# Description: 用于在完成全部流程后，给文档重命名，之后或许可以在源代码中增添这一步，从而去除当前脚本
# ==== 主体代码开始 ====
while read -r old_val new_val; do
    for file in *"${old_val}"*.pdf; do
        if [ -f "$file" ]; then
            new_name="${file/${old_val}/sim${new_val}}"
            mv "$file" "$new_name"
            echo "Renamed: $file -> $new_name"
        fi
    done
done < rename.txt

#rename.txt
# sim1	0.001
# sim2	0.049
# sim3	0.349
# sim4	0.586
# sim5	0.664
# sim6	0.8
# sim7	3.691

# Renamed: QT_Amy_2_scatter_sim1.pdf -> QT_Amy_2_scatter_sim0.001.pdf
