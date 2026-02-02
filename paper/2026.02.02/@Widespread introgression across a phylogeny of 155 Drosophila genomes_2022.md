---
title: Widespread introgression across a phylogeny of 155 Drosophila genomes
authors: Anton Suvorov, Bernard Y. Kim, Jeremy Wang, Ellie E. Armstrong, David Peede, Emmanuel R. R. D’Agostino, Donald K. Price, Peter J. Waddell, Michael Lang, Virginie Courtier-Orgogozo, Jean R. David, Dmitri Petrov, Daniel R. Matute, Daniel R. Schrider, Aaron A. Comeault
year: 2022
citekey: suvorovWidespreadIntrogressionPhylogeny2022
tags:
  - paper
  - literature
  - 果蝇
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Widespread introgression across a phylogeny of 155 Drosophila genomes
</div>

**Authors:** Anton Suvorov, Bernard Y. Kim, Jeremy Wang, Ellie E. Armstrong, David Peede, Emmanuel R. R. D’Agostino, Donald K. Price, Peter J. Waddell, Michael Lang, Virginie Courtier-Orgogozo, Jean R. David, Dmitri Petrov, Daniel R. Matute, Daniel R. Schrider, Aaron A. Comeault  
**Year:** 2022  
**Zotero:** [Open in Zotero](zotero://select/items/@suvorovWidespreadIntrogressionPhylogeny2022)

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

</div>

---
##  文献使用的数据集

155个果蝇的基因组数据
使用BUSCOs筛选获得了2791个单拷贝直系同源数据集

## 文献使用的主要方法
系统发育上：
- IQTREE与ASTRAL
分化时间估算上：
- MCMCtree
- BEAST2
杂交渐渗上;
- DCT **Discordant Count Test** 不一致基因计数测试      自设计方法
- BLT **Branch Length Test** 枝长测试      自设计方法
- QuIBL
- Dsuite
- HyDe
- ==**f-branch 启发式方法**== 可以用于解缠和映射在多个相关物种对中检测到的渗入事件。它通过整合DCT/BLT估算的渗入比例 (γ)，来估算整个演化支系中渗入事件的数量和基因组的渗入比例
- **PhyloNet**



## 文献的主要内容与理论
1） 核心问题：
- 果蝇属中存在的杂交渐渗情况
- 杂交渐渗主要发生在哪些内群中
- 杂交渐渗强度如何
- 发生的时间是古代还是近代


2）主要结论：
- **广泛存在的基因渗入** ，**DCT和BLT的证据**: 在分析的9个演化支中，有7个演化支至少有一个物种对同时被DCT和BLT检测出显著的渗入证据（FDR校正p值<0.05）。在演化支1和3中，尽管两个测试并非同时显著，但都指向相同的渗入物种对
- **古老和近期渗入并存**: 将基因树不一致性映射到系统发育树上显示，在所检查的9个演化支中的大多数都发生了古老和近期的渗入，揭示了基因流在整个果蝇生命树中的持续作用。
- **渗入基因组比例**，f-branch启发式方法估计了不同基因组的基因渗入比例。无论渗入事件的发生时间早晚，物种间交换的基因组比例（γ）相对相似，这表明渗入并非总是有害且被快速清除。
- **PhyloNet**分析，与前面的杂交渐渗分析得到了几乎一致的结论，使得渗入检测结果更加可靠和保守

3）理论创新点：
1.**祖先种群结构（ancestral population structure）的某些情景可能会导致基因树在两种不一致拓扑结构下的数量和支长出现差异** 。在物种分化之前，祖先种群内部可能存在的分层或亚群结构。如果祖先种群存在亚群，并且这些亚群之间的基因流不均匀，那么即使没有发生物种间的渗入，仅仅是祖先种群内部的结构，也可能在后代物种的基因组中留下类似渗入的信号。

解决方案：祖先种群结构主要影响基因树的合并时间深度，可能导致不一致拓扑之间的合并时间差异，从而模拟渗入信号。但祖先种群结构通常不会导致不一致拓扑的合并时间比一致拓扑的合并时间更近期（浅）。如果一个基因流事件是真的渗入（发生在物种分化之后），它应该导致不一致拓扑的合并时间比物种树上同源基因的合并时间（一致拓扑的支长）更近期。通过与一致拓扑的支长进行比较，可以更有效地过滤掉那些只是祖先种群结构造成的深度合并时间差异，而不能反映近期基因流的信号。
**比较了不一致拓扑基因树的支长与一致拓扑基因树的支长**

2.根据鉴定出的最有可能发生渐渗的物种对，选定进行phylonetwork分析的类群。并能够验证此前的假设


## 思考

往后的研究可以结合
- 谱系地理
- 全基因
- 泛基因组
来进行展开

这样的优势在于，现在同一个物种的种群中，其基因型可能还是存在一定的差异，这就使得取样上的差异可能会导致最终显示出基因树的不一致，而这种不一致可能与杂交渐渗无关，而是反映了种群不同个体基因型的差异。因此从种群的角度和泛基因组的角度考虑可能能够古更有效精准的解决这一问题。

全基因组的优势不言而喻。基因取样上更加具有代表性，不存在说由于取样的单拷贝直系同源基因更多的集中在杂交区域而导致过高的估计渐渗基因在基因组中的占比。

