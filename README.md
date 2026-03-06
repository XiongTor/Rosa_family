# 蔷薇科系统发育相关研究
# Rosaceae Phylogenomics Analysis Pipeline

本项目主要包含研究中所涉及的代码脚本，涵盖数据获取、数据过滤、系统发育分析、多样化分析、祖先重建及网状进化与ILS水平估计等分析流程。

This repository contains comprehensive scripts for phylogenomic analysis of the Rosaceae family, including data acquisition, phylogenetic tree construction, incongruence analysis (ILS, GTEE, IH), introgression detection, and distribution data analysis.

**作者 (Author)**: Tao Xiong

---

## 研究流程图 (Research Workflow) - 2025.02

![蔷薇科果实系统发育研究](https://github.com/XiongTor/Rosa_family/blob/main/imag/readme_flowwork/%E8%94%B7%E8%96%87%E7%A7%91%E6%9E%9C%E5%AE%9E%E7%B1%BB%E5%9E%8B%E7%B3%BB%E7%BB%9F%E5%8F%91%E8%82%B2%E7%A0%94%E7%A9%B6.svg)

![蔷薇科多倍化及网状进化研究](https://github.com/XiongTor/Rosa_family/blob/main/imag/readme_flowwork/%E8%94%B7%E8%96%87%E7%A7%91%E7%BD%91%E7%8A%B6%E8%BF%9B%E5%8C%96%E5%8F%8A%E4%B8%8D%E5%AE%8C%E5%85%A8%E8%B0%B1%E7%B3%BB%E7%AD%9B%E9%80%89%E7%A0%94%E7%A9%B6.svg)

![蔷薇科核型演化研究](https://github.com/XiongTor/Rosa_family/blob/main/imag/readme_flowwork/%E8%94%B7%E8%96%87%E7%A7%91%E6%A0%B8%E5%9E%8B%E6%BC%94%E5%8C%96.svg)

---

## 目录 (Table of Contents)

- [数据下载 (Data Download)](#数据下载-data-download)
  - [性状数据 (Character Data)](#性状数据-character-data)
  - [分布数据 (Distribution Data)](#分布数据-distribution-data)
- [分析 (Analysis)](#分析-analysis)
  - [建树分析 (Tree Building)](#建树分析-tree-building)
  - [Hybpiper分析 (Hybpiper Analysis)](#hybpiper分析-hybpiper-analysis)
  - [ILS分析 (ILS Analysis)](#ils分析-ils-analysis)
  - [BAMM分析 (BAMM Analysis)](#bamm分析-bamm-analysis)
- [工具 (Tools)](#工具-tools)
- [统计分析与绘图 (Statistical Analysis)](#统计分析与绘图-statistical-analysis)
- [其他 (Others)](#其他-others)

---

## 数据下载 (Data Download)

### 性状数据 (Character Data)
用于从各种数据库下载和处理性状数据的脚本。

| 脚本 (Script) | 描述 (Description) | 作者 (Author) |
|--------|-------------|--------|
| [get_dispersal_type.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/character_data_down/get_dispersal_type.R) | 基于扩散术语对每个物种的扩散类型进行分类，并探索扩散类型与分布范围之间的关系 <br> Classify dispersal type and explore relationship with distribution range | Tao Xiong |
| [get_character_data.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/character_data_down/get_character_data.R) | 从TRY、BIEN、GIFT数据库下载性状数据并合并 <br> Download character data from TRY, BIEN, GIFT databases | Tao Xiong |

### 分布数据 (Distribution Data)
用于下载、清理和可视化物种分布数据的脚本。

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [get_distribution_data.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/distribution_data_down/get_distribution_data.R) | 从GBIF、idigbio、CVH、NSII数据库获取蔷薇科分布数据 <br> Get distribution data from GBIF, idigbio, CVH, NSII | Tao Xiong | 2025.03.14 |
| [get_wcvp_native_distribution.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/distribution_data_down/get_wcvp_native_distribution.R) | 合并数据集，标准化物种名称，通过WCVP获取原生分布 <br> Combine datasets, standardize names, get native distribution via WCVP | Tao Xiong | 2025.03.14 |
| [clean_distribution_data.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/distribution_data_down/clean_distribution_data.R) ⚙️ | 清理分布数据（需要Species、Longitude、Latitude列）<br> Clean distribution data (requires Species, Longitude, Latitude columns) | Tao Xiong | 2025-02-26 |
| [map_rosa_subf.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/distribution_data_down/map_rosa_subf.R) | 将蔷薇亚科坐标信息映射到世界地图 <br> Map Rosa subfamily coordinates to world map | Tao Xiong | 2025.03.15 |
| [count_continent_num_by_coordinate.R](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/distribution_data_down/count_continent_num_by_coordinate.R) | 计算各大陆原生和引入物种的数量 <br> Calculate native and introduced species per continent | Tao Xiong | 2025.03.15 |
| [change_longtitude_to_address_baidu.py](https://github.com/XiongTor/Rosa_family/blob/main/code/DATA_DOWNLOAD/distribution_data_down/change_longtitude_to_address_baidu.py) | 使用百度API将经纬度转换为地址 <br> Reverse geocoding using Baidu API | Tao Xiong | 2025.03.14 |

> ⚙️ 标记表示该脚本可直接在命令行运行 (Executable scripts)

---

## 分析 (Analysis)

### 建树分析 (Tree Building)
使用多种方法（RAxML、ASTRAL、IQ-TREE等）进行系统发育树构建的脚本。

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [Pyphlawd_analysis.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/Pyphlawd_analysis.sh) | 设置PyPHLAWD数据库并为蔷薇科系统发育分析寻找优质聚类 <br> Setup PyPHLAWD database for Rosaceae phylogenetic analysis | Tao Xiong | 2025.03.27 |
| [make_tree.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/make_tree.sh) | 多种建树方法记录（modeltest、RAxML-ng等）<br> Various tree building methods (modeltest, RAxML-ng) | - | - |
| [genetree_reroot.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/genetree_reroot.sh) | 使用外类群或MAD方法重定根基因树 <br> Reroot gene trees using outgroup or MAD method | Tao Xiong | 2025-08-27 |
| [pyphlawd_semi_auto.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/pyphlawd_semi_auto.sh) | 使用PyPHLAWD进行半自动系统发育树构建 <br> Semi-automatic phylogenetic tree building using PyPHLAWD | Tao Xiong | 2025-12-16 |
| [hybsuite_use.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/hybsuite_use.sh) | 使用hybsuite（目标捕获流程）构建系统发育树 <br> Build tree using hybsuite (target capture pipeline) | Tao Xiong | 2025-08-27 |
| [treeshrink_pipeline.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/treeshrink_pipeline.sh) | 过滤长度小于100bp的基因序列并进行treeshrink分析 <br> Filter sequences <100bp and run treeshrink analysis | Tao Xiong | 2026-01-11 |
| [plot_tree.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/plot_tree.R) | 蔷薇科细胞核项目的自动化树绘图 <br> Auto plot trees for Rosaceae cytonuclear project | Tao Xiong | 2026-02-09 |
| [rosa_make_tree_old_2022.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/bulidtree/rosa_make_tree_old_2022.sh) | 2022年的旧建树流程（存档）<br> Old tree-building pipeline from 2022 (archive) | xiongtao | 2023.05.16 |

### Hybpiper分析 (Hybpiper Analysis)
用于目标捕获数据组装和树构建的脚本。

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [trimmomatic.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/hybpiper_analy/trimmomatic.sh) | 使用Trimmomatic对Illumina测序数据进行质量修剪 <br> Quality trimming using Trimmomatic | Tao Xiong | 2025-03.26 |
| [hybpiper_step.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/hybpiper_analy/hybpiper_step.sh) | 目标捕获数据的逐步hybpiper组装流程 <br> Step-by-step hybpiper assembly process | Tao Xiong | 2025-03.26 |
| [make_tree_step.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/hybpiper_analy/make_tree_step.sh) | 记录hybpiper分析建树步骤（MAFFT比对、trimAL修剪）<br> Record tree building steps (MAFFT, trimAL) | Tao Xiong | 2025-04.04 |
| [hybpiper_paftol_combine.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/hybpiper_analy/hybpiper_paftol_combine.sh) | 合并hybpiper和PAFTOL 353基因数据集 <br> Combine hybpiper and PAFTOL 353-gene datasets | - | - |

### ILS分析 (ILS Analysis)
用于分析不完全谱系分选（ILS）、基因树估计误差（GTEE）和渐渗杂交（IH）的综合脚本。

#### 主要ILS脚本 (Main ILS Scripts)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [QS_test.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/QS_test.sh) | ILS分析的四重奏采样（QS）测试 <br> Quartet Sampling (QS) test for ILS analysis | Tao Xiong | 2025-09-26 |
| [phytop.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/phytop.sh) | 运行ASTRAL和phytop进行物种树可视化 <br> Run ASTRAL and phytop for species tree visualization | Tao Xiong | 2025-10-17 |
| [count_sp_rich_value.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/count_sp_rich_value.sh) | 过滤后计算每个基因的物种丰富度值 <br> Count species richness per gene after filtering | Tao Xiong | 2026-02-09 |
| [phybase_simulate.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/phybase_simulate.R) | 使用phybase模拟不同的进化场景 <br> Simulate evolutionary scenarios using phybase | Tao Xiong | 2025-10-24 |

#### 细胞核显示 (Cytonuclear Display)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [QS_visual.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/cytonuclear_display/QS_visual.R) | 可视化四重奏采样（QS）结果 <br> Visualize Quartet Sampling results | Xiongtao | - |
| [RF_path_MDS.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/cytonuclear_display/RF_path_MDS.R) | Robinson-Foulds距离的MDS分析 <br> MDS analysis of Robinson-Foulds distances | Tao Xiong | 2025-10-18 |
| [densitree.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/cytonuclear_display/densitree.R) | 绘制densitree以显示基因树之间的差异 <br> Plot densitree to show gene tree distinction | Tao Xiong | 2025-07-22 |
| [phyparts.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/cytonuclear_display/phyparts.sh) | 运行phyparts分析并生成饼图可视化 <br> Run phyparts and generate pie chart visualization | - | - |
| [phyparts_pies.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/cytonuclear_display/phyparts_pies.sh) | 将phyparts分析结果显示为饼图 <br> Display phyparts results as pie charts | - | - |

#### D统计分析 (D-statistics Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [D-statistics.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/D-statistics.sh) | 计算D统计量以检查基因渐渗 <br> Calculate D-statistics to check gene introgression | Tao Xiong | 2025-05-16 |
| [run_Dsuite.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/run_Dsuite.py) | Dsuite Dtrios自动化流程（枚举群体、提取序列、运行分析）<br> Auto pipeline for Dsuite Dtrios | Tao Xiong | 2025-10-14 |
| [run_Dsuite_xt.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/run_Dsuite_xt.py) | 四群体测试的Dsuite Dfoil流程 <br> Dsuite Dfoil pipeline for 4-population test | Tao Xiong | 2025-10-14 |
| [run_Dsutie_xt_snp_thinning.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/run_Dsutie_xt_snp_thinning.py) | 带SNP稀疏化的Dsuite流程 <br> Dsuite pipeline with SNP thinning | Tao Xiong | 2025-10-14 |
| [run_Dsuite_no_snp_thinning.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/run_Dsuite_no_snp_thinning.py) | 不带SNP稀疏化的Dsuite流程 <br> Dsuite pipeline without SNP thinning | Tao Xiong | 2025-10-14 |
| [msa2vcf_biallelic.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/msa2vcf_biallelic.py) | 从MSA fasta提取双等位基因SNP并输出VCF <br> Extract biallelic SNPs from MSA and output VCF | Tao Xiong | 2026-01-16 |
| [run.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/run.sh) | 使用nohup并行运行所有生成的Dsuite脚本 <br> Run all Dsuite scripts in parallel with nohup | Tao Xiong | 2026-01-16 |
| [make_runbash.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/make_runbash.sh) | 为每个参数组合生成单独的bash脚本 <br> Generate individual bash scripts for each parameter combination | Tao Xiong | 2026-01-16 |
| [plot_box.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/D-statistics/plot_box.R) | 亚科间和亚科内D/f4值的比较 <br> Comparison of D/f4 values between/within subfamilies | Tao Xiong | 2025-10-16 |

#### ILS/GTEE/IH分析 (ILS/GTEE/IH Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [get_node_inf.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/ILS_GTEE_IH/get_node_inf.R) | 从树中获取节点信息 <br> Get node information from trees | Tao Xiong | 2025-07-04 |
| [get_theta_value.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/ILS_GTEE_IH/get_theta_value.R) | 从ASTRAL树和重新缩放的分支长度树获取theta值 <br> Get theta values from ASTRAL and rescaled trees | Tao Xiong | 2025-07-04 |
| [relaimpo.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/ILS_GTEE_IH/relaimpo.R) | 计算ILS、GTEE和IH因子的相对重要性 <br> Calculate relative importance of ILS, GTEE, IH | Tao Xiong | 2025-07-04 |

#### MSCquartets分析 (MSCquartets Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [MSCquartets.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/MSCquartets/MSCquartets.R) | 记录MSCquartets分析工作流程 <br> Record MSCquartets analysis workflow | Tao Xiong | 2025.05.04 |
| [get_gene_ILS_IH_contribution.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/MSCquartets/get_gene_ILS_IH_contribution.R) | 评估ILS和IH对基因树不一致性的贡献 <br> Evaluate ILS and IH contributions to gene tree incongruence | Tao Xiong | 2026.01.28 |
| [get_gene_ILS_IH_contribution_old.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/MSCquartets/get_gene_ILS_IH_contribution_old.R) | 评估ILS和IH的贡献（旧版本）<br> Evaluate ILS and IH contributions (older version) | Tao Xiong | 2025.11.14 |
| [get_sp_ILS_IH_contribution.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/MSCquartets/get_sp_ILS_IH_contribution.R) | 获取物种水平的ILS和IH贡献 <br> Get species-level ILS and IH contributions | - | - |
| [results_classify.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/MSCquartets/results_classify.sh) | 对MSCquartets分析结果进行分类 <br> Classify results from MSCquartets analysis | Tao Xiong | 2026-01-29 |
| [analy_df_num.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/MSCquartets/analy_df_num.R) | 分析results_classify.sh的结果（冲突分析和可视化）<br> Analyze results from results_classify.sh | Tao Xiong | 2026-02-06 |

#### 系统发育网络分析 (PhyloNetwork Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [PhyloNetworks.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/PhyloNetwork/PhyloNetworks.R) | 运行PhyloNetworks（SnaQ）分析进行系统发育网络推断 <br> Run PhyloNetworks (SnaQ) for network inference | Tao Xiong | 2025-05-15 |
| [plot_result_julia.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/PhyloNetwork/plot_result_julia.R) | 显示PhyloNetworks分析结果（基于Julia的可视化）<br> Display PhyloNetworks results (Julia-based) | Tao Xiong | 2026-02-12 |

#### QuIBL分析 (QuIBL Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [run_quibl.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/QuIBL/run_quibl.sh) | 运行QuIBL分析 <br> Run QuIBL analysis | Tao Xiong | 2026-01-16 |
| [QuIBL_plot.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/QuIBL/QuIBL_plot.R) | 绘制QuIBL分析结果 <br> Plot QuIBL analysis results | Tao Xiong | 2026-02-08 |

#### 模拟研究 (Simulation Studies)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [phybase_simulate.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/simulate/phybase_simulate.R) | 使用phybase模拟不同的进化场景 <br> Simulate evolutionary scenarios using phybase | Tao Xiong | 2025-10-24 |
| [get_theta.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/simulate/get_theta.R) | 通过MRCA匹配边来计算theta，然后注释树 <br> Compute theta by matching edges via MRCA | Tao Xiong | 2025-10-30 |
| [count_four_aray.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/simulate/count_four_aray.R) | 从基因树中计算四重奏频率 <br> Count quartet frequencies from gene trees | Tao Xiong | 2025.05.04 |
| [count_four_aray_lm.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/simulate/count_four_aray_lm.R) | 从基因树中计算四重奏频率（线性模型版本）<br> Count quartet frequencies (linear model version) | Tao Xiong | 2025.11.18 |

#### 因子比对树分析 (Factor Alignment Tree Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [factor.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/factor_alignment_tree/factor.sh) | 使用phykit和AMAS获取比对和树的统计信息 <br> Use phykit and AMAS to get alignment/tree statistics | Tao Xiong | 2025-12-07 |
| [factor_parallel.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/factor_alignment_tree/factor_parallel.sh) | 因子比对树分析的并行处理 <br> Parallel processing for factor alignment tree analysis | Tao Xiong | - |
| [make_tree_factor.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/factor_alignment_tree/make_tree_factor.sh) | 为因子比对的20次迭代构建树 <br> Build trees for 20 iterations of factor alignment | Tao Xiong | 2025-12-16 |
| [factor_plot.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/factor_alignment_tree/factor_plot.R) | 分析蔷薇科叶绿体和核基因组之间的不同因子 <br> Analyze factors between Chloroplast and Nuclear genomes | XiongTao | 2025.12.11 |

#### 树外类群添加 (Tree Outgroup Addition)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [add_outgroup_gene_specific_change.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/tree_addoutg/add_outgroup_gene_specific_change.py) | 根据最近的内群亲缘关系添加具有基因特异性分支长度的外类群 <br> Add outgroup with gene-specific branch lengths | Tao Xiong | 2026-01-16 |
| [keep_one_outg.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/tree_addoutg/keep_one_outg.sh) | 在基因树中仅保留一个外类群 <br> Keep only one outgroup in gene trees | Tao Xiong | 2026-01-16 |

#### IH基因获取 (Get ILS IH Gene)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [get_and_plot_IH_gene.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/ILS/Get_ILS_IH_gene/get_and_plot_IH_gene.R) | 读取IH/ILS基因结果并绘制箱线图和条形图 <br> Read IH/ILS gene results and plot boxplots/barplots | Tao Xiong | 2025-11-19 |

### BAMM分析 (BAMM Analysis)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [plot_bamm_env.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Analysis/BAMM/plot_bamm_env.R) | BAMM相关绘图 <br> BAMM related plotting | Tao Xiong | 2025-05-13 |

---

## 工具 (Tools)

用于各种数据处理和分析任务的实用脚本。

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [wcvp.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/wcvp.R) ⚙️ | 通过WCVP进行名称标准化 <br> Name standardization by WCVP | Tao Xiong | 2025.03.28 |
| [mono.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/mono.R) ⚙️ | 检查每个分支的单系性 <br> Check monophyly of each clade | Tao Xiong | 2025.04.16 |
| [monophy.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/monophy.R) ⚙️ | 使用R包MonoPhy检查单系性 <br> Check monophyly using MonoPhy package | Tao Xiong | 2025.04.18 |
| [make_hybpiper_reference.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/make_hybpiper_reference.sh) | 获取353基因数据的参考序列 <br> Get reference sequence for 353-gene data | Tao Xiong | 2025-04.03 |
| [genetree_reroot.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/genetree_reroot.sh) | 使用外类群或MAD方法重定根基因树 <br> Reroot gene trees using outgroup or MAD | Tao Xiong | 2025-08-27 |
| [count_OG_per_species.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/count_OG_per_species.sh) | 计算每个物种的直系同源组（OGs）数量 <br> Count orthologous groups per species | Tao Xiong | 2025-07-02 |
| [clean_supp.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/clean_supp.R) | 清理支持值 <br> Clean support values | Tao Xiong | 2025-07-03 |
| [count_RF_path_dist.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/count_RF_path_dist.R) | 计算Robinson-Foulds和路径距离 <br> Count Robinson-Foulds and path distances | - | - |
| [count_RF_path_two_trees.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/count_RF_path_two_trees.R) | 计算两棵树之间的RF和路径距离 <br> Count RF and path distance between two trees | Tao Xiong | 2025-10-08 |
| [msa2vcf_biallelic.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/msa2vcf_biallelic.py) | 将MSA FASTA转换为双等位基因VCF格式 <br> Convert MSA FASTA to biallelic VCF | - | - |
| [gcf_scf.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/gcf_scf.sh) | 计算基因和位点一致性因子 <br> Count gene and site concordance factors | Tao Xiong | 2025-10-23 |
| [remove_lt_100bp_length_moregap.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/remove_lt_100bp_length_moregap.sh) | 移除短于100bp的序列 <br> Remove sequences shorter than 100bp | Tao Xiong | 2025-10-11 |
| [make_treelength_gt_0.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/make_treelength_gt_0.R) ⚙️ | 使树的分支长度大于0 <br> Make tree branch lengths greater than 0 | Tao Xiong | 2025.04.18 |
| [split_windows.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/split_windows.sh) | 将序列分割为滑动窗口 <br> Split sequences into sliding windows | Tao Xiong | 2025-12-02 |
| [pair_distance_aa_dna.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/pair_distance_aa_dna.R) | 序列距离分析 <br> Sequence distance analysis | Tao Xiong | 2025-12-21 |
| [kaks.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/kaks.sh) | 使用pal2nal计算Ka/Ks比率 <br> Calculate Ka/Ks ratios using pal2nal | Tao Xiong | 2025-12-21 |
| [count_length_gap.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/count_length_gap.sh) | 计算序列长度和gap统计信息 <br> Count sequence length and gap statistics | - | - |
| [Get_support_value.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/Get_support_value.R) | 提取带有节点ID的ASTRAL支持值 <br> Extract ASTRAL support values with node IDs | Tao Xiong | 2025-08-27 |
| [cophylo.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/cophylo.R) | Cophylo分析 - 比较物种树的拓扑结构 <br> Cophylo analysis - compare species tree topology | Tao Xiong | 2025-12-18 |
| [GDR.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/GDR.py) | 使用Selenium从rosaceae.org JBrowse抓取数据 <br> Scrape data from rosaceae.org JBrowse | - | - |
| [pdf_to_png.py](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/pdf_to_png.py) | 将PDF转换为具有目标宽度的PNG <br> Convert PDF to PNG with target width | Tao Xiong | 2026-03-04 |
| [reboot_frps_smoothy.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Tools/reboot_frps_smoothy.sh) | 重启FRPC以进行远程服务器连接 <br> Reboot FRPC for remote server connection | - | - |

> ⚙️ 标记表示该脚本可直接在命令行运行 (Executable scripts)

**使用示例 (Usage Examples):**

```bash
# WCVP名称标准化
Rscript wcvp.R name.csv

# 检查单系性
Rscript mono.R test.tree
Rscript monophy.R test.tre

# 使树长度大于0
Rscript make_treelength_gt_0.R test.tre

# 清理分布数据
Rscript clean_distribution_data.R data/distribution_data.csv
```

---

## 统计分析与绘图 (Statistical Analysis and Graphing)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [Boxplot.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Statistical_Analysis_and_Graphing/Boxplot.R) | 绘制箱线图（模板/示例）<br> Plot boxplot (template/example) | Tao Xiong | 2025-12-18 |

---

## 其他 (Others)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) |
|--------|-------------|--------|
| [Zijia_build_tree.sh](https://github.com/XiongTor/Rosa_family/blob/main/code/Others/Zijia_build_tree.sh) | 从353基因中解析序列并按属组织 <br> Parse sequences from 353 genes and organize by genus | Miao Sun |
| [poaceae_distributon_wcvp_plot.R](https://github.com/XiongTor/Rosa_family/blob/main/code/Others/poaceae_distributon_wcvp_plot.R) | 绘制黑麦草和其他黑麦草属物种的分布 <br> Plot distribution of Lolium perenne and other Lolium | Tao Xiong |

---

## 图片处理 (Picture Processing)

| 脚本 (Script) | 描述 (Description) | 作者 (Author) | 日期 (Date) |
|--------|-------------|--------|------|
| [pdf_to_png.py](https://github.com/XiongTor/Rosa_family/blob/main/code/picture/pdf_to_png.py) | 将PDF转换为PNG图像 <br> Convert PDF to PNG images | Tao Xiong | 2026-02-02 |
| [combine_png.py](https://github.com/XiongTor/Rosa_family/blob/main/code/picture/combine_png.py) | 将文件夹中的所有PDF合并为4x6网格图像 <br> Merge all PDFs into a 4x6 grid image | Tao Xiong | 2026-02-02 |

---

## 统计摘要 (Summary Statistics)

- **Shell脚本总数 (Total Shell Scripts)**: 37
- **R脚本总数 (Total R Scripts)**: 45
- **Python脚本总数 (Total Python Scripts)**: 13
- **脚本总数 (Total Scripts)**: 95

---

## 主要分析重点 (Key Analysis Focus Areas)

1. **系统发育树构建 (Phylogenetic Tree Building)**: RAxML, ASTRAL, IQ-TREE, MPEST等多种方法
2. **不一致性分析 (Incongruence Analysis)**: ILS、GTEE、IH测量
3. **渐渗检测 (Introgression Detection)**: D统计量、Dsuite
4. **四重奏分析 (Quartet Analysis)**: 四重奏采样、MSCquartets
5. **分布分析 (Distribution Analysis)**: GBIF、WCVP、物种分布制图
6. **序列处理 (Sequence Processing)**: 比对、修剪、质量控制
7. **树比较 (Tree Comparisons)**: Robinson-Foulds距离、路径距离、拓扑比较
8. **网络分析 (Network Analysis)**: PhyloNetworks、QuIBL

---

## 引用 (Citation)

如果您使用这些脚本，请引用分析中使用的相应工具和方法。

If you use these scripts, please cite the appropriate tools and methods used in your analysis.

---

## 联系方式 (Contact)

**Tao Xiong**
GitHub: [@XiongTor](https://github.com/XiongTor)
Email: 1760755795@qq.com

---

*最后更新 (Last updated): 2026-03-05*