# RF Distance Analysis - Two Strategies

本目录包含两种不同的等距区间划分策略用于RF距离分析。

## 策略对比

### 策略 4.1：基于密度分布的等距区间（Density-based）
**脚本**: `08_density_based_rf_analysis.R`

**原理**：
- 根据属性值的**分布密度**来划分区间
- 在数值空间上选取相等长度的区间
- 每个区间覆盖相同的数值范围，但包含的基因数量可能不同

**特点**：
- 区间边界由数值范围决定（如 0.40-0.45, 0.50-0.55）
- 适合分析数值分布特征
- 不同区间的基因数量可能差异较大
- 更关注数值空间的均匀性

**示例**：
```r
# 设置分位点（成对出现，每对定义一个区间）
quantile_points <- c(0.40, 0.45, 0.5, 0.55, 0.605, 0.65, 0.71, 0.76, 0.82, 0.86)
# 结果：每个区间覆盖相同的数值范围比例
```

---

### 策略 4.2：基于数量的等距区间（Quantile-based）
**脚本**: `05_rf_distance_analysis.R`

**原理**：
- 根据基因**数量**来划分区间
- 将基因按属性值排序后，选取相同数量的基因作为一个区间
- 每个区间包含相同数量的基因，但数值范围可能不同

**特点**：
- 区间边界由基因排序位置决定（如前10%-15%的基因）
- 适合比较不同区间的统计特性
- 每个区间的基因数量相等或接近
- 更关注样本数量的均衡性

**示例**：
```r
# 设置分位点
quantile_points <- c(0.1, 0.15, 0.25, 0.3, 0.4, 0.45, 0.55, 0.6, 0.7, 0.75)
# 结果：Low组包含10%-15%的基因，Low_m组包含25%-30%的基因
```

---

## 使用方法

### 运行策略 4.1（密度分布）
```bash
cd E:\Project\Rosacea\code\Analysis\ILS\factor_alignment_tree\factor_plot
Rscript 08_density_based_rf_analysis.R
```

**配置参数**：
```r
factor <- "Mean_internal_branch"  # 修改为要分析的因子
interval_mode <- "paired"          # "paired" 或 "consecutive"
quantile_points <- c(0.40, 0.45, 0.5, 0.55, 0.605, 0.65, 0.71, 0.76, 0.82, 0.86)
```

### 运行策略 4.2（数量等分）
```bash
Rscript 05_rf_distance_analysis.R
```

**配置参数**：
```r
factor <- "Mean_internal_branch"  # 修改为要分析的因子
interval_mode <- "paired"          # "paired" 或 "consecutive"
quantile_points <- c(0.1, 0.15, 0.25, 0.3, 0.4, 0.45, 0.55, 0.6, 0.7, 0.75)
```

---

## 输出文件

### 策略 4.1 输出
- `{factor}_density_based_equidistants.pdf`
  - 密度分布图（带高亮区间）
  - 因子值分布的小提琴图
  - RF距离分布的小提琴图

### 策略 4.2 输出
- `{factor}_equidistants.pdf`
  - 基因排序柱状图（带分组标记）
  - 因子值分布的小提琴图
  - RF距离分布的小提琴图

---

## 选择建议

**使用策略 4.1（密度分布）当**：
- 需要研究特定数值范围的基因特征
- 关注数值分布的形态特征
- 想要在数值空间上均匀采样

**使用策略 4.2（数量等分）当**：
- 需要比较不同组的统计特性
- 希望每组有相似的样本量以提高统计检验力
- 关注排序位置而非绝对数值

---

## 依赖关系

两个脚本都依赖：
- `00-utils_functions.R` - 工具函数（包含 `pairwise_rf_onecol`）
- `merged_orthofinder.rds` - 预处理的数据
- `./tree/orthofinder_genes/` - 基因树文件

确保先运行 `01_data_loading.R` 生成必要的数据文件。
