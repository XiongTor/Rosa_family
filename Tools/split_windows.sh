#!/bin/python3
# Author: Tao Xiong
# Date: 2025-12-02
# Description: split windows
# ==== 主体代码开始 ====
#!/bin/bash

# 1) 滑窗（200bp，步长100，最小长度100）
seqkit sliding -W 300 -s 150 -g cPrSOPT4_CDS.trim.fasta -o tmp_sliding.fasta

# 2) 重组窗口
grep "^>" tmp_sliding.fasta | sed 's/.*://' | sort -u | while read win; do
    out="window_${win}.fasta"
    awk -v w="$win" '
        /^>/ { keep = ($0 ~ ":"w"$") }
        keep { print }
    ' tmp_sliding.fasta > "$out"
done
