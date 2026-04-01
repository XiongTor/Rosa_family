---
title: "An integrative framework reveals widespread gene flow during the early radiation of oaks and relatives in Quercoideae (Fagaceae)"
authors: Shui‐Yin Liu, Ying‐Ying Yang, Qin Tian, Zhi‐Yun Yang, Shu‐Feng Li, Paul J. Valdes, Alex Farnsworth, Heather R. Kates, Carolina M. Siniscalchi, Robert P. Guralnick, Douglas E. Soltis, Pamela S. Soltis, Gregory W. Stull, Ryan A. Folk, Ting‐Shuang Yi
year: 2025
citekey: liuIntegrativeFrameworkReveals2025
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  An integrative framework reveals widespread gene flow during the early radiation of oaks and relatives in Quercoideae (Fagaceae)
</div>

**Authors:** Shui‐Yin Liu, Ying‐Ying Yang, Qin Tian, Zhi‐Yun Yang, Shu‐Feng Li, Paul J. Valdes, Alex Farnsworth, Heather R. Kates, Carolina M. Siniscalchi, Robert P. Guralnick, Douglas E. Soltis, Pamela S. Soltis, Gregory W. Stull, Ryan A. Folk, Ting‐Shuang Yi  
**Year:** 2025  
**Zotero:** [Open in Zotero](zotero://select/items/@liuIntegrativeFrameworkReveals2025)

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
  Although the frequency of ancient hybridization across the Tree of Life is greater than previously thought, little work has been devoted to uncovering the extent, timeline, and geographic and ecological context of ancient hybridization. Using an expansive new dataset of nuclear and chloroplast DNA sequences, we conducted a multifaceted phylogenomic investigation to identify ancient reticulation in the early evolution of oaks (Quercus). We document extensive nuclear gene tree and cytonuclear discordance among major lineages of Quercus and relatives in Quercoideae. Our analyses recovered clear signatures of gene flow against a backdrop of rampant incomplete lineage sorting, with gene flow most prevalent among major lineages of Quercus and relatives in Quercoideae during their initial radiation, dated to the Early‐Middle Eocene. Ancestral reconstructions including fossils suggest ancestors of Castanea + Castanopsis, Lithocarpus, and the Old World oak clade probably co‐occurred in North America and Eurasia, while the ancestors of Chrysolepis, Notholithocarpus, and the New World oak clade co‐occurred in North America, offering ample opportunity for hybridization in each region. Our study shows that hybridization—perhaps in the form of ancient syngameons like those seen today —has been a common and important process throughout the evolutionary history of oaks and their relatives. Concomitantly, this study provides a methodological framework for detecting ancient hybridization in other groups.
</div>

---

## 1. 文献研究类群与使用的数据集

#### **研究类群：**
稻属
11属27种，包含6个二倍体属（AA, BB, CC, EE, FF, GG）和5个多倍体属（BBCC, CCDD, HHJJ,  HHKK, KKLL）

#### **原始数据集：**
**自测数据**包含**基因组重测序数据**与**转录组数据**
同时抓取了二倍体属的网上公共数据集的**染色体级别基因组数据**。

#### **直系同源基因数据集：**
通过OrthoFinder筛选组装出的包含**3973 single-copy orthologs**

#### **SNP数据获取**


## 2. 文献使用的主要方法
#### 建树方法：
串联法建树：RAxML建树，分别使用SNP和单拷贝直系同源基因建树
溯祖法建树：RAxML建立基因树，ASTRALⅢ建立物种树

同时也建立了质体树：使用RAxML建立

#### 基因树不一致分析方法： 
Phyparts
SplitsTree
PhyloNet
D-statistics (Dfm与Dxy): Dfm滑动窗口确定具体基因位置
QuIBL： BIC显著性判断ILS，IH占比
GO富集
TreeMix（推断历史群体层面基因流动的方向和规模）
r8s （推断分化时间）




## 3. 文献的主要内容与理论
**1） 核心问题**：稻属早期研究难以得到一致的令人信服的系统发育结果。前人有提出稻属存在严重的杂交渐渗现象，但是并没有具体量化。  
**2）主要结论：**   
- 稻属内存在强烈的杂交渐渗现象，在杂交渐渗和ILS中，文章的结论支持杂交渐渗占据主导地位，约占总比的74.17%。  
- 通过染色体级别基因组数据，鉴定出具体的渐渗基因，通过GO富集发现这类基因主要富集在环境适应（应激反应、次级代谢）相关的功能上，证明了网状进化在水稻适应环境中的关键作用
 
**3） 理论创新点：** 解决了如何量化杂交渐渗的问题，同时鉴定出了具体的渐渗基因，有助于了解杂交渐渗在物种适应与演化中起到的作用。同时也为后续类似研究提供了范例
   

## 4. 思考
文章中目前最大的问题是，似乎忽略了**系统误差带来的干扰**，对于IH占比74.17的结论，几乎是只通过QuIBL的显著率来判断，而没有具体分析其ILS的可能情况，可能存在ILS强度低但是显著的情况，而这种低水平的ILS可能就是背景噪音造成的。同理IH的高占比也没有核实并作出相应解释。




