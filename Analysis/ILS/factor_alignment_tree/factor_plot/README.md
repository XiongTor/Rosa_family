# Rosaceae Cytonuclear Factor Analysis

## 脚本说明

原始脚本 `factor_plot.R` 已被拆分为多个模块化脚本，便于单独运行和维护。

## 文件结构

```
├── utils_functions.R          # 所有自定义函数
├── 01_data_loading.R          # 数据读取和预处理
├── 02_density_plots.R         # 密度图绘制
├── 03_correlation_analysis.R  # 相关性分析
├── 04_pca_analysis.R          # PCA分析
├── 05_rf_distance_analysis.R  # RF距离分析
├── 06_chloroplast_matching.R  # 叶绿体相似基因匹配
├── 07_final_comparison.R      # 最终比较分析
└── run_all.R                  # 主控脚本
```

## 使用方法

### 方式1: 运行完整流程（自动化）

```r
source("run_all.R")
```

这将自动运行步骤1-4，步骤5-7需要手动配置。

### 方式2: 单独运行各个模块

#### 步骤1: 数据加载
```r
source("01_data_loading.R")
```
- 读取叶绿体和核基因数据
- 生成 `merged_chloroplast.rds` 和 `merged_orthofinder.rds`

#### 步骤2: 密度图
```r
source("02_density_plots.R")
```
- 生成所有因子的密度分布图
- 输出: `all_density_plots.pdf`

#### 步骤3: 相关性分析
```r
source("03_correlation_analysis.R")
```
- 计算因子间相关性
- 输出: `corrplot_diff_all_factor.pdf`, `col.pdf`

#### 步骤4: PCA分析
```r
source("04_pca_analysis.R")
```
- 执行主成分分析
- 输出: `pca.pdf`, `df_all.rds`, `pca.rds`

#### 步骤5: RF距离分析（需配置）
```r
# 编辑脚本中的参数:
# - factor: 要分析的因子名称
# - quantile_points: 分位点设置
source("05_rf_distance_analysis.R")
```
- 计算不同因子区间的RF距离
- 输出: `{factor}_equidistants.pdf`

#### 步骤6: 叶绿体相似基因匹配（需配置）
```r
# 编辑脚本中的参数:
# - factor: 要匹配的因子
# - pct_value: 匹配容差百分比
source("06_chloroplast_matching.R")
```
- 找出与叶绿体基因相似的核基因
- 输出: CSV文件和密度图

#### 步骤7: 最终比较（需要树文件）
```r
source("07_final_comparison.R")
```
- 比较不同因子的RF距离
- 输出: `Boxplot_all_diff.pdf`

## 主要改进

1. **模块化**: 每个分析步骤独立，便于调试和重复运行
2. **数据持久化**: 使用RDS文件保存中间结果，避免重复计算
3. **代码简化**: 移除冗余代码，保持核心功能
4. **灵活配置**: 关键参数集中在脚本开头，便于修改

## 依赖包

```r
library(ape)
library(phangorn)
library(dplyr)
library(tidyverse)
library(Hmisc)
library(corrplot)
library(cowplot)
library(ggplot2)
library(plotly)
```

## 注意事项

- 步骤5-7需要根据具体分析需求手动配置参数
- 确保数据文件路径正确（chloroplast/, orthofinder/, tree/）
- 步骤7需要先构建系统发育树
