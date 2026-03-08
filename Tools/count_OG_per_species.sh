#!/bin/bash
# Author: Tao Xiong
# Date: 2025-07-02
# Description: This script is used to count the number of orthologous groups (OGs) per species from Orthogroups/Orthogroups.GeneCount.tsv
# ==== 主体代码开始 ====
#!/bin/bash

# 输入文件

input_file=$1
name=$(echo $1|cut -f2 -d'/')
output_file=${name}_species_OG_num.csv

nn=$(ls *.csv |grep -c $name)
if [ $nn -gt 0 ];then
    echo "$output_file is exit"
    num=$((nn + 1))
    output_file=${name}_${num}_species_OG_num.csv
fi

# 检查输入文件是否存在
if [ ! -f "$input_file" ]; then
    echo "Error: $input_file not found!"
    exit 1
fi

# 获取表头（物种名）
header=$(head -n 1 "$input_file")
IFS=$'\t' read -r -a species <<< "$header"

# 初始化输出
echo -e "Species,OG_Count,single_OG" > "$output_file"

# 对每个物种列（从第2列开始）统计非零行数
num_species=${#species[@]}
for (( i=1; i<num_species; i++ )); do
    sp=${species[$i]}
    count=$(awk -F '\t' -v col=$((i+1)) 'NR>1 && $col > 0 {c++} END{print c+0}' "$input_file")
    single=$(awk -F '\t' -v col=$((i+1)) 'NR>1 && $col == 1 {c++} END{print c+0}' "$input_file")
    echo -e "${sp},${count},${single}" >> "$output_file"
done

echo "Done. Result saved to $output_file"