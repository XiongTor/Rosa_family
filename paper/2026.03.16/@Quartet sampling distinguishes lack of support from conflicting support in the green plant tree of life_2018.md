---
title: "Quartet sampling distinguishes lack of support from conflicting support in the green plant tree of life"
authors: James B. Pease, Joseph W. Brown, Joseph F. Walker, Cody E. Hinchliff, Stephen A. Smith
year: 2018
citekey: peaseQuartetSamplingDistinguishes2018
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Quartet sampling distinguishes lack of support from conflicting support in the green plant tree of life
</div>

**Authors:** James B. Pease, Joseph W. Brown, Joseph F. Walker, Cody E. Hinchliff, Stephen A. Smith  
**Year:** 2018  
**Zotero:** [Open in Zotero](zotero://select/items/@peaseQuartetSamplingDistinguishes2018)

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
  Premise of the Study Phylogenetic support has been difficult to evaluate within the green plant tree of life partly due to a lack of specificity between conflicted versus poorly informed branches. As data sets continue to expand in both breadth and depth, new support measures are needed that are more efficient and informative. Methods We describe the Quartet Sampling ( QS ) method, a quartet‐based evaluation system that synthesizes several phylogenetic and genomic analytical approaches. QS characterizes discordance in large‐sparse and genome‐wide data sets, overcoming issues of alignment sparsity and distinguishing strong conflict from weak support. We tested QS with simulations and recent plant phylogenies inferred from variously sized data sets. Key Results QS scores demonstrated convergence with increasing replicates and were not strongly affected by branch depth. Patterns of QS support from different phylogenies led to a coherent understanding of ancestral branches defining key disagreements, including the relationships of Ginkgo to cycads, magnoliids to monocots and eudicots, and mosses to liverworts. The relationships of ANA ‐grade angiosperms ( Amborella , Nymphaeales, Austrobaileyales), major monocot groups, bryophytes, and fern families are likely highly discordant in their evolutionary histories, rather than poorly informed. QS can also detect discordance due to introgression in phylogenomic data. Conclusions Quartet Sampling is an efficient synthesis of phylogenetic tests that offers more comprehensive and specific information on branch support than conventional measures. The QS method corroborates growing evidence that phylogenomic investigations that incorporate discordance testing are warranted when reconstructing complex evolutionary histories, in particular those surrounding ANA ‐grade, monocots, and nonvascular plants.
</div>

---

## 1. 核心问题
#### 1）系统发育支持度的评估困难，难以区分部分冲突位点是信息不足还是进化历史本身冲突
#### 2）针对部分"**大型稀疏矩阵**"数据集，现有的主流UFboot和PP值难以准确评估其系统发育支持度
UFboot或者说传统的 Bootstrap（NBS）在对于”大型稀疏矩阵“（gap过多的矩阵）时，容易多次抽取到无信息的位点
PP即贝叶斯法在面对大型数据时计算速率过慢



## 2. 主要方法
Quartet sampling 基于==**四联体**==的方法来处理上述核心问题
**核心机理**：更具选取的位点所连接的内部分支，将系统发育树拆分成四个主要的部分
评估某一节点的QC/QD/QI值时，是从四个部分中随机挑选一个来共同组成一个quartet。然后针对这个quartet去评估它所有三种可能的拓扑中似然值最高的一个拓扑为当前运算的最优拓扑，如果此最优quartet所展现出的拓扑结构与系统发育树相吻合，则为这个节点记入一次支持。
- ### QC
     (Quartet Concordance，四分体一致性)：取值 [-1, 1]。衡量抽样四分体中有多少支持主拓扑。1表示主拓扑占据绝对优势
- ### QD
	  (Quartet Differential，四分体差异/冲突不对称性)：取值 [0, 1]。如果不支持主树的四分体，在另外两种备选拓扑中分布**不均匀**（偏斜，QD值低），通常暗示存在**基因流/杂交/渐渗**；如果分布均匀（QD值接近1），则更可能是 ILS 造成的随机噪音。
- ### QI
     (Quartet Informativeness，四分体信息量)：取值 [0, 1]。评估这些四分体是否具有足够的系统发育信号（区分拓扑的能力），即最优拓扑的对数似然值显著高于次优拓扑的对数似然值，一般默认差值在2以上即有显著差异
![](../imag/@Quartet%20sampling%20distinguishes%20lack%20of%20support%20from%20conflicting%20support%20in%20the%20green%20plant%20tree%20of%20life_2018/file-20260312221206513.png)



# 3. 对于杂交渐渗的鉴定作用
通过模拟数据证实了：
- ILS升高，会导致QC降低，QD升高
- IH升高，QC，QD均降低

# 4. 为系统发育节点的稳固性和不一致性提供新视角

1） 对于数据量丰富的数据集而言，其BS值大概率维持在高值，但是QC值可能为负，暗示了其物种树拓扑存在争议，可能存在复杂的演化历史

**ANA 基部类群（ANA-grade）** 的演化关系
利用 QS 对四个主流数据集（WI2014, ZN2014, HS2014, XI2014）进行了检测，发现所有数据集中关于次分支的QC值都为负值，说明可能所有的研究都无法给出十分确凿的证据

2）**对于不同数据集而言，虽然它们的物种树指向了不同的结果，但QC值可能暗示了相同的进化历史**;

**真双子叶植物（Eudicots）**、**单子叶植物（Monocots）** 和 **木兰类（Magnoliids）** 之间的演化先后顺序
**主流假说 ：** 支持“木兰类 + 真双子叶植物”组成姐妹群，而单子叶植物较早分化
**非主流假说：** 显示“真双子叶植物 + 单子叶植物”是姐妹群。
 QS 发现非主流假说使用的数据集在这个分支上的 QC 值为负。这意味着虽然物种树画出了该分支，但实际上从quartet出发可能更多地投向了另一种可能。
 经过 QS 评估，原本矛盾的四个数据集达成了一致——它们其实都支持**木兰类与真双子叶植物具有更近的共同祖先**，而排除了单子叶植物