#!/bin/bash
# Author: Tao Xiong
# Date: 2026-01-29
# Description: Classify the results obtained before analysis.
# ==== 主体代码开始 ====


# 1.数据分类 -----------------------------------------------------------------------------------

# 按照大的分类分出不同的文件夹，并将对应的文件移入相应的文件夹中
mkdir df_num pairs_summary conflict_matrix plot_pdf
mv *df_num* df_num
mv *pairs_summary* pairs_summary
mv *conflict_matrix* conflict_matrix
mv *pdf plot_pdf/

# 如果文件太长可以使用：
find . -maxdepth 1 -name "*pairs_summary*" -print0 \
| xargs -0 mv -t pairs_summary



# 按照不同的alpha与beta组合来划分文件夹
# 对df_num文件夹内的文件不需要进行分类，其它的文件夹下可以进行一定的分类
alpha_list=("0.01" "0.001" "1e-04" "1e-05" "1e-06")
beta_list=("0.01" "0.05")
data_list=("red" "yellow" "green")

# 定义alpha的标准化名称（用于文件夹命名）
declare -A alpha_names
alpha_names["0.01"]="0.01"
alpha_names["0.001"]="0.001"
alpha_names["1e-04"]="0.0001"
alpha_names["1e-05"]="0.00001"
alpha_names["1e-06"]="0.000001"

# 按照alpha、beta和data组合创建子文件夹并移动文件
for alpha in "${alpha_list[@]}"; do
    for beta in "${beta_list[@]}"; do
        for data in "${data_list[@]}"; do
            # 使用标准化的alpha名称创建文件夹
            folder_name="${data}_alpha${alpha_names[$alpha]}_beta${beta}"
            
            # 检查是否有匹配的文件
            if ls *${data}*alpha_${alpha}_beta_${beta}* 2>/dev/null | grep -q .; then
                mkdir -p "$folder_name"
                mv *${data}*alpha_${alpha}_beta_${beta}* "$folder_name"/
                echo "已移动 ${data} alpha=${alpha} beta=${beta} 的文件到 ${folder_name}/"
            fi
        done
    done
done

# 2. df_num文件夹文件合并与分析 --------------------------------------------------------------------
mkdir -p orgin_data & mv * orgin_data/

for file in orgin_data/*.csv;do
    cat $file|sed '1d'>>all_df_num.csv
done

sed -i "1i gene,all_quartet,all_diff,blue,red,yellow,green,alpha_value,beta_value" all_df_num.csv

# 若处理的是orthofinder的数据，需要把gene名从 OG0013459_1_noout_output_rt_oneoutg_final 改为 OG0013459_1
awk -F',' 'BEGIN{OFS=","} NR==1{print; next}
{
  split($1, a, "_")
  $1 = a[1] "_" a[2]
  print
}' all_df_num.csv > all_df_num_gene2.csv

sed -i "s/\"//g" all_df_num_gene2.csv