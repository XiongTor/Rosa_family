---
title: "Demographic history and natural selection shape patterns of deleterious mutation load and barriers to introgression across <i>Populus</i> genome"
authors: Shuyu Liu, Lei Zhang, Yupeng Sang, Qiang Lai, Xinxin Zhang, Changfu Jia, Zhiqin Long, Jiali Wu, Tao Ma, Kangshan Mao, Nathaniel R Street, Pär K Ingvarsson, Jianquan Liu, Jing Wang
year: 2022
citekey: liuDemographicHistoryNatural2022
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Demographic history and natural selection shape patterns of deleterious mutation load and barriers to introgression across <i>Populus</i> genome
</div>

**Authors:** Shuyu Liu, Lei Zhang, Yupeng Sang, Qiang Lai, Xinxin Zhang, Changfu Jia, Zhiqin Long, Jiali Wu, Tao Ma, Kangshan Mao, Nathaniel R Street, Pär K Ingvarsson, Jianquan Liu, Jing Wang  
**Year:** 2022  
**Zotero:** [Open in Zotero](zotero://select/items/@liuDemographicHistoryNatural2022)

---

<div style="
  border-radius: 8px;
  background-color: #2f2d2d;
  padding: 15px 20px;
  margin-top: 15px;
  color: #ddd6d6;
  line-height: 1.6;
  font-size: 16px;
">
  <div style="font-size: 22px;font-weight: bold; color: #bbb3b3; margin-bottom: 8px;">
    ❝ Abstract
  </div>
  Hybridization and resulting introgression are important processes shaping the tree of life and appear to be far more common than previously thought. However, how the genome evolution was shaped by various genetic and evolutionary forces after hybridization remains unresolved. Here we used whole-genome resequencing data of 227 individuals from multiple widespread Populus species to characterize their contemporary patterns of hybridization and to quantify genomic signatures of past introgression. We observe a high frequency of contemporary hybridization and confirm that multiple previously ambiguous species are in fact F1 hybrids. Seven species were identified, which experienced different demographic histories that resulted in strikingly varied efficacy of selection and burdens of deleterious mutations. Frequent past introgression has been found to be a pervasive feature throughout the speciation of these Populus species. The retained introgressed regions, more generally, tend to contain reduced genetic load and to be located in regions of high recombination. We also find that in pairs of species with substantial differences in effective population size, introgressed regions are inferred to have undergone selective sweeps at greater than expected frequencies in the species with lower effective population size, suggesting that introgression likely have higher potential to provide beneficial variation for species with small populations. Our results, therefore, illustrate that demography and recombination have interplayed with both positive and negative selection in determining the genomic evolution after hybridization.
</div>

---

## 1. 文献研究类群与使用的数据集
杨属(*Populus*)
采样227个个体，获得全基因组重测序数据


## 2. 文献使用的主要方法
#### Call SNP
选取*P. tremula*的全基因组数据为参考，其余数据使用32×重测序，进行Call  SNP，获得了8,966,513个SNP
使用NJ法建树



## 3. 文献的主要内容与理论
#### 核心问题：
杂交发生后，进入受体物种的“外来基因片段”面临怎样的命运？是什么进化力量（如遗传漂变、自然选择、重组率）决定了哪些片段被保留，哪些被淘汰？
#### 文章目的：
利用227个个体，12种的杨属植物个体全基因组重测序数据，表征杨属的当代杂交模式，并量化过去渗入的基因组特征。
#### 整体路线
（1）先用全基因组数据理清物种分类和当代杂交状态；
（2）比较不同物种的种群历史和遗传负荷差异；
（3）定位基因组上的历史渗入区域；
（4）最后，将渗入区域与重组率、遗传负荷、选择信号进行关联分析，揭示选择力量的作用机制。

#### （1）先用全基因组数据理清物种分类和当代杂交状态
通过SNP数据进行PCA聚类分析，表明了有7个物种可以明显的被分离开，另外5个物种被认为是这7个物种的杂交后代

通过祖源推断与fastSTRUCTURE，叶绿体，线粒体系统发育树以及杂合性和杂交指数表明：许多过去被认为是一个独立“种”的杨树，实际上是不同物种间自发产生的**杂交种**。

| **杂交种名称                    | **母本 (Maternal)**       | **父本 (Paternal)**       | **杂交程度                  |
| -------------------------- | ----------------------- | ----------------------- | ----------------------- |
| **毛白杨** *`P. tomentosa`*   | 响叶杨 *`P. adenopoda`*    | 银白杨 *`P. alba`*         | 自发产生的 $F_1$ 代杂交种        |
| **西藏杨** *`P. tibetica`*    | 圆叶杨 *`P. rotundifolia`* | 银白杨 *`P. alba`*         | $F_1$ 代杂交种              |
| **五莲杨** *`P. wulianensis`* | 山杨 *`P. davidiana`*     | 响叶杨 *`P. adenopoda`*    | $F_1$ 代杂交种              |
| **宁陕杨** *`P. ningshanica`* | 山杨 *`P. davidiana`*     | 响叶杨 *`P. adenopoda`*    | $F_1$ 代杂交种（母本源自不同的山杨群体） |
| **银灰杨** *`P. canescens`*   | 欧洲山杨 *`P. tremula`*     | 银白杨 *`P. alba`*         | 包含 $F_1$ 代及与亲本回交的后代     |
| **山杨×圆叶杨杂交区**              | 山杨 *`P. davidiana`*     | 圆叶杨 *`P. rotundifolia`* | 存在大量具有混合基因型的杂交个体        |

#### （2) 评估和比较七个物种的遗传多样性和种群历史
移除了上述5个混和种后，使用剩下的7个物种的162个个体进行群体遗传分析
综合评估了:
- 系统发育结果
- **遗传分化**$F_{ST}$ 与 $d_{xy}$   ：衡量物种间的遗传距离和分化程度
- **遗传多样性**核苷酸多样性 ($\pi$)、单体型 SNP (Singletons) ：评估种群内的遗传变异丰富度
- **连锁不平衡**LD 衰减分析 ($r^2$ 衰减速度) ：反映重组率和种群历史大小（LD 衰减越慢，通常有效种群越小）
- **种群历史推断** **PSMC** (Pairwise Sequential Markovian Coalescent)

得到如下结论：

| **物种 / 类别**                  | **关键发现**                                 | **演化结论**                                   |
| ---------------------------- | ---------------------------------------- | ------------------------------------------ |
| **琼岛杨** *`P. qiongdaoensis`* | 多样性极低，单体型 SNP 最少，LD 衰减最慢。                | **长期受困**：由于岛屿隔离，种群规模长期严重萎缩，面临较高的灭绝风险或遗传瓶颈。 |
| **响叶杨** *`P. adenopoda`*     | 多样性较低，PSMC 显示在末次冰盛期（LGM）后有所增长。           | **先降后升**：经历过长时间衰退，但近期环境转好后种群开始复苏。          |
| **美洲颤杨** *`P. tremuloides`*  | 单体型 SNP 极多，Tajima's D 为显著负值，PSMC 显示急剧扩张。 | **成功扩张**：在冰期结束后经历了剧烈的种群增长，占据了广泛的生境。        |
| **银白杨** *`P. alba`*          | PSMC 轨迹独特，约 0.05-0.07 Ma 达到峰值后持续下降。      | **独特轨迹**：其生存历史与其他杨树不同，可能受特定的地理或气候事件影响。     |
| **高多样性组** (山杨等)              | 保持了较高水平的核苷酸多样性 ($\pi$)。                  | **稳健发展**：这些物种分布广泛，受环境波动影响相对较小，保持了丰富的遗传资源。  |
![](../imag/@Demographic%20history%20and%20natural%20selection%20shape%20patterns%20of%20deleterious%20mutation%20load%20and%20barriers%20to%20introgression%20across%20_i_Populus__i_%20genome_2022/file-20260420152632416.png)

#### 评估了有害突变负荷
主要通过以下几点来评估：

|**维度**|**方法 / 工具**|**目的**|
|---|---|---|
|**选择压力估计**|**DFE-alpha** 软件|估算不同强度的突变（从中性到强有害）在物种中的比例。|
|**近交与繁殖模式**|**$F_{IS}$ (近交系数)**|评估物种是否存在近交或克隆生长。|
|**突变负荷量化**|**0/4倍简并位点杂合度比率**|通过比较非同义与同义变异的比例，衡量清除有害突变的效率。|
|**突变功能预测**|**SIFT4G** 软件|将变异分类为：同义、可耐受、有害、功能缺失 (LoF)。|
|**等位基因极化**|使用近缘种作为外类群|确定哪个碱基是祖先态，哪个是新产生的（衍生态）。|

发现琼岛杨的近中性突变 ($N_e s < 1$) 比例显著升高并表现出极负的 $F_{IS}$ 和高杂合度。说明该种群太小，主要靠**克隆繁殖**（无性系）维持。能够保护了杂合度，但进一步阻碍了有害突变的重组清除。

**在小群体中，由于选择压力松弛，这些致命变异躲在杂合子中不表现出表型，从而逃避了自然选择的清除，被长期保留下来**


#### （3）定位基因组上的历史渗入区域------==主要关注内容==
1. 通过TWISST来评判全基因组中不同的滑动窗口的拓扑占比。发现占比最高的拓扑与此前的NJ树的结果一致，最后进行了ABBA分析，通过D值和f4值来观察各个物种间的杂交渐渗情况。
2. 考虑到琼岛杨可能的克隆繁殖现象，因此后续不将其纳入考量范围内
![|675](../imag/@Demographic%20history%20and%20natural%20selection%20shape%20patterns%20of%20deleterious%20mutation%20load%20and%20barriers%20to%20introgression%20across%20_i_Populus__i_%20genome_2022/file-20260420163430965.png)

3. 随后将剩下的6个物种，组成10个具有渐渗现象的trio，通过滑动窗口计算他们各自的Fdm值，其中滑动窗口大小为50个SNPs，步长为20SNPs。然后确定Fdm值的阈值，即高于多少才算是渐渗区段。首先计算全局的渐渗比例，例如为x%.则取Fdm值全X%的区域为渐渗区段。之后从多个指标评判了渐渗区段：
- 物种间分化度和物种内分化度(DXY，FST)，发现渐渗区段拥有更低的物种间分化，受体物种内拥有更高的物种内多样性。
- 重组率，发现渗入区域显著集中在高重组区域，且其包含的有害突变与同义变异的比例比随机水平低
- 基因功能，远缘物种与生殖功能相关，近缘物种与应激反应，DNA修复相关


##### 适应性渐渗，判断哪些区段是适应性的渐渗区段：
我们使用VolcanoFinder计算了复合似然比（CLR）统计量，以检测每个物种内的渗入扫描特征(signatures of introgression sweeps)。适应性内渗区域定义为与推断中最高5% CLR支持的内渗扫荡重叠的区域。即这部分区域收到了强烈的自然选择。发现：
- 在有效群体大小较小的类群中，适应性渗入显著，而在有效群体大小较大的类群中，适应性渗入不显著
- 高编码区密度（基因很多）和低重组率区域可能会导致适应性渗入的假阳性，而检测到的适应性渗入区段多包含非编码区，说明不存在假阳性的可能
- 通过单倍型检测 iHH12 发现，适应性渗入区域的 iHH12 值显著更高，再次验证了其真实性。  iHH12：一种基于单倍型纯合度（Haplotype Homozygosity）的检测方法。如果一个有利变异在群体中迅速扩散，它周围的单倍型会非常整齐、单一。
- GO富集，发下适应性渐渗区段的功能主要集中在生殖调节，生物钟，抗逆反应，分子信号上






