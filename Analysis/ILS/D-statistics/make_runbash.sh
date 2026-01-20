#!/bin/bash
# Author: Tao Xiong
# Date: 2026-01-16
# Description: Generate individual bash scripts for each parameter combination

# 定义参数数组
min_cov_values=(0.04 0.5 0.6 0.7 0.8 0.9)
min_maf_values=(0.05 0.01)

echo "Generating scripts for all parameter combinations..."
echo "Total combinations: $((${#min_cov_values[@]} * ${#min_maf_values[@]}))"
echo "----------------------------------------"

# 遍历所有参数组合
for cov in "${min_cov_values[@]}"; do
    for maf in "${min_maf_values[@]}"; do
        # 创建对应的文件夹和脚本名称
        folder_name="Dsuite_${cov}_${maf}"
        script_name="Dsuite_${cov}_${maf}.sh"

        # 创建文件夹
        mkdir -p "${folder_name}"

        # 生成bash脚本内容（直接放在文件夹内）
        cat > "${folder_name}/${script_name}" << EOF
#!/bin/bash
# Author: Tao Xiong
# Date: 2026-01-16
# Description: D-statistics calculation with min_cov=${cov}, min_maf=${maf}

# ==== 主体代码开始 ====

# 获取vcf文件
python ~/data/scripts/msa2vcf_biallelic.py ../rosa_MO_orthofinder_supermatrix.fasta rosa_orthofinder.vcf --min_cov ${cov} --min_maf ${maf}

########################## D-statistics calculation ##############################

Dsuite Dtrios rosa_orthofinder.vcf ../sets.txt -t ../rosa_orthofinder_MO_treeshrink_sp_rt.tre -o result

ruby ~/data/scripts/Dsuite/plot_d.rb result_tree.txt ../plot_order.txt 0.7 species_sets_no_geneflow_BBAA_D.svg

ruby ~/data/scripts/Dsuite/plot_f4ratio.rb result_tree.txt ../plot_order.txt 0.2 species_sets_no_geneflow_BBAA_f4ratio.svg

# 计算f-branch值

Dsuite Fbranch ../rosa_orthofinder_MO_treeshrink_sp_rt.tre result_tree.txt >fbranch.out


# 用dtools.py脚本绘制f-branch图

dtools.py fbranch.out ../rosa_orthofinder_MO_treeshrink_sp_rt.tre --outgroup Zelkova_schneideriana --use_distances --dpi 500 --tree-label-size 20

echo "Completed: min_cov=${cov}, min_maf=${maf}"
EOF

        # 添加执行权限
        chmod +x "${folder_name}/${script_name}"

        echo "✓ Generated: ${folder_name}/${script_name}"
    done
done

echo "----------------------------------------"
echo "All scripts generated successfully!"
echo ""
echo "To run all scripts in parallel, use:"
echo "  for folder in Dsuite_*; do (cd \"\$folder\" && bash *.sh) & done; wait"
echo ""
echo "Or run them one by one:"
echo "  cd Dsuite_0.04_0.01 && bash Dsuite_0.04_0.01.sh"
