---
title: "Phylogenomic analyses reveal species relationships and phylogenetic incongruence with new member detected in allium subgenus cyathophora"
authors: Kun Chen, Zi-Jun Tang, Yuan Wang, Jin-Bo Tan, Song-Dong Zhou, Xing-Jin He, Deng-Feng Xie
year: 2025
citekey: chenPhylogenomicAnalysesReveal2025
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Phylogenomic analyses reveal species relationships and phylogenetic incongruence with new member detected in allium subgenus cyathophora
</div>

**Authors:** Kun Chen, Zi-Jun Tang, Yuan Wang, Jin-Bo Tan, Song-Dong Zhou, Xing-Jin He, Deng-Feng Xie  
**Year:** 2025  
**Zotero:** [Open in Zotero](zotero://select/items/@chenPhylogenomicAnalysesReveal2025)

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
  Species characterized by undetermined clade affiliations, limited research coverage, and deficient systematic investigation serve as enigmatic entities in plant and animal taxonomy, yet hold critical significance for exploring phylogenetic relationships and evolutionary trajectories. Subgenus Cyathophora (Allium, Amayllidaceae), a small taxon comprising approximately five species distributed in the Qinghai–Tibet Plateau (QTP) and adjacent regions might contain an enigmatic species that has long remained unexplored. In this study, we collected data on species from subgenus Cyathophora and its close relatives in subgenus Rhizirideum, as well as the enigmatic species Allium siphonanthum. Combining phylogenomic datasets and morphological evidence, we investigated species relationships and the underlying mechanism of phylogenetic discordance. A total of 1662 single-copy genes (SCGs) and 150 plastid loci were filtered and used for phylogenetic analyses based on concatenated and coalescent-based methods. Furthermore, to systematically evaluate phylogenetic discordance and decipher its underlying drivers, we implemented integrative analyses using multiple approaches, such as coalescent simulation, Quartet Sampling (QS), and MSCquartets. Our phylogenetic analyses robustly resolve A. siphonanthum as a member of subg. Cyathophora, forming a sister clade with A. spicatum. This relationship was further corroborated by their shared morphological characteristics. Despite the robust phylogenies inferred, extensive phylogenetic conflicts were detected not only among gene trees but also between SCGs and plastid-derived species trees. These significant phylogenetic incongruences in subg. Cyathophora predominantly stem from incomplete lineage sorting (ILS) and reticulate evolutionary processes, with historical hybridization events likely correlated with the past orogenic dynamics and paleoclimatic oscillations in the QTP and adjacent regions. Our findings not only provide new insights into the phylogeny of subg. Cyathophora but also significantly enhance our understanding of the evolution of species in this subgenus.
</div>

---

## 1. 文献研究类群与使用的数据集
- 研究类群为**葱属杯花亚属（*Allium* subgenus *Cyathophora*）**
- 收集了12种葱属植物，进行转录组测序
- orthofinder产生1662个单拷贝基因
- GetOrganelle产生150个质体基因位点数据（包括质体基因间隔区，IGSs）


## 2. 文献使用的主要方法
- IQTREE---串联
- ASTRALⅢ---溯祖
- Quartet Sampling
- MSCquartets
- Phytop
- PhyloNet
- NANUQ
- **基因一致性因子（gCF）和位点一致性因子（sCF）**



## 3. 文献的主要内容与结论
1） 核心问题：论文旨在解决葱属（Allium L.）杯花亚属（subgenus Cyathophora）中物种关系不明确的问题，特别是针对一个长期被忽视、系统学位置存疑的物种**多花葱（*Allium siphonanthum*）**。并探究**系统发育不一致性及其驱动因素**
2）主要结论：
-  **多花葱（Allium siphonanthum）被确认为杯花亚属的新成员**，并与A. spicatum形成姐妹群，这解决了其长期未知的分类学位置。
- **不完全谱系分选（ILS）是造成系统发育不一致性的主要原因**，其信号在多个关键节点上表现强烈。
- **网状演化（杂交）是次要但重要的驱动因素**，研究识别出A. tetraploideum和A. siphonanthum均存在杂交起源，表明杂交在亚属进化中扮演了关键角色。
- **历史地质构造运动和古气候振荡**（如青藏高原的隆升和第四纪气候波动）可能促进了这些ILS和杂交过程，从而推动了物种分化和环境适应。


## 4. 主要关注点--- gcf，scf，QS的数据如何分析

1）单独分析核基因和质体基因的gcf，scf大小
- **核基因组数据（SCGs）**：gCF 值（38.5%~95.4%，平均 65.89%）**始终高于** sCF 值（38.0%~88.9%，平均 59.24%）。这表明在细胞核基因组中，大多数完整的基因树是支持物种树拓扑的。
- **质体基因组数据**：情况相反。gCF 值极低（16.7%~64.7%，平均 34.96%），但 sCF 值很高（42.4%~95.1%，平均 71.08%）。（注：这在生物学上是合理的，因为质体基因组通常作为一个整体不重组遗传，强行将其拆分成单个基因建树，由于许多短基因缺乏足够的变异位点，会导致单基因树拓扑混乱，从而拉低 gCF；而看所有变异位点 sCF 时，一致性又很高。）

2）计算gcf，scf值与BS，枝长的相关性
- 总体而言呈现正相关，即枝长越长，gcf，scf值越大。BS值始终较为稳定
- 出现了**短分支，低gcf,scf，高BS**的情况，暗示了ILS的可能性（出现了多个拓扑不同但支持率较高的情况）

3）关于QS检测
- 主要是通过定位**高BS，高QC，低QD**的节点，明确它们暗示了潜在的杂交网络
- 通过QF值来明确取样的准确性，不存在“流氓类群”
