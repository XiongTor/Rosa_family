---
title: "Novel information theory-based measures for quantifying incongruence among phylogenetic trees"
authors: Leonidas Salichos, Alexandros Stamatakis, Antonis Rokas
year: 2014
citekey: salichosNovelInformationTheorybased2014
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Novel information theory-based measures for quantifying incongruence among phylogenetic trees
</div>

**Authors:** Leonidas Salichos, Alexandros Stamatakis, Antonis Rokas  
**Year:** 2014  
**Zotero:** [Open in Zotero](zotero://select/items/@salichosNovelInformationTheorybased2014)

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
  Phylogenies inferred from different data matrices often conflict with each other necessitating the development of measures that quantify this incongruence. Here, we introduce novel measures that use information theory to quantify the degree of conflict or incongruence among all nontrivial bipartitions present in a set of trees. The first measure, internode certainty (IC), calculates the degree of certainty for a given internode by considering the frequency of the bipartition defined by the internode (internal branch) in a given set of trees jointly with that of the most prevalent conflicting bipartition in the same tree set. The second measure, IC All (ICA), calculates the degree of certainty for a given internode by considering the frequency of the bipartition defined by the internode in a given set of trees in conjunction with that of all conflicting bipartitions in the same underlying tree set. Finally, the tree certainty (TC) and TC All (TCA) measures are the sum of IC and ICA values across all internodes of a phylogeny, respectively. IC, ICA, TC, and TCA can be calculated from different types of data that contain nontrivial bipartitions, including from bootstrap replicate trees to gene trees or individual characters. Given a set of phylogenetic trees, the IC and ICA values of a given internode reflect its specific degree of incongruence, and the TC and TCA values describe the global degree of incongruence between trees in the set. All four measures are implemented and freely available in version 8.0.0 and subsequent versions of the widely used program RAxML.
</div>

---

## 1. 核心问题

传统的BS指标以及MRC（majority-rule consensus）树存在一定缺陷，无法得知与主拓扑相互冲突的拓扑是怎样的分布频率

# 2. 核心方法
## 2.1 二分法
### **定义**：一个内部节点（分支）会将所有物种切分成两个互不相交的集合，这就叫二分法。
同时每一个包含有K个物种的二歧树都将包含k-3个**非平凡二分 (nontrivial bipartitions)**
非平凡二分即：即这 $k-3$ 个二分中的每一个都将树中的 $k = m + n$ 个类群划分为两个分别包含 $m$ 和 $n$ 个类群的分组，其中 $m \ge 2$ 且 $n \ge 2$

### 二分的**相容与冲突**
来自同一类群集合（Taxon set）的两个二分 $A = X_1 | X_2$ 和 $B = Y_1 | Y_2$ 被称为“**相容的**”，当且仅当这四个子集对的交集中至少有一个为空集：$(X_1 \cap Y_1, X_1 \cap Y_2, X_2 \cap Y_1, X_2 \cap Y_2)$成立，如果不成立，则认为它们是**相互冲突**的

让我们考虑二分 $A = \{a, b, c, d, e \mid f, g, h, i, j\}$，它由子集 $A_1 = \{a, b, c, d, e\}$ 和 $A_2 = \{f, g, h, i, j\}$ 组成（其中 $a$ 到 $j$ 是类群名称）。

再考虑同一类群集中的第二个二分 $B = \{a, b, c \mid d, e, f, g, h, i, j\}$，它由子集 $B_1 = \{a, b, c\}$ 和 $B_2 = \{d, e, f, g, h, i, j\}$ 组成。**二分 $B$ 与二分 $A$ 不冲突**，因为 $A_2 \cap B_1$ 是空集。

相比之下，第三个二分 $C = \{a, b, c, d, g \mid e, f, h, i, j\}$，由子集 $C_1 = \{a, b, c, d, g\}$ 和 $C_2 = \{e, f, h, i, j\}$ 组成。**二分 $C$ 与二分 $A$ 冲突（或不相容）**，因为四个交集 $(A_1 \cap C_1, A_1 \cap C_2, A_2 \cap C_1, A_2 \cap C_2)$ 均不为空

# 2.2 信息论与香农熵

### 香农熵：是衡量一个系统中**信息量大小**或**不确定性程度**的指标，熵越大越乱，随机性越大，信息量越大
香农熵的计算公式如下：
$$H(X) = -\sum_{i=1}^{n} P(x_i) \log_b P(x_i)$$
举一个例子：
当我们抛掷硬币的时候，正反面出现的概率都为0.5，此时香农熵最大，不确定性最大，结果最难以预测。
而当硬币被做了手脚，正面出现的概率为0.99的时候，此时香农熵会减小，不确定性很小，结果很容易预测。

文中香农熵的概念引入到了系统发育树分支冲突的评估中来，并以此判断物种树种某一的节点拓扑的确定性
IC只考虑主拓扑与主次要拓扑
ICA考虑主拓扑与其它所有拓扑，但从计算的角度，一般包含占比在替代拓扑中5%以上的替代拓扑
其中：IC 和 ICA的计算公式和函数图如下：
$$IC = 1 + P(X_1) \log_2[P(X_1)] + P(X_2) \log_2[P(X_2)]$$  

<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.03.16/file-20260314161308147.png width="60%">


$$ICA = \log_n(n) + \sum_{i=1}^{n} P(X_i) \log_n[P(X_i)]$$

<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.03.16/file-20260314161457213.png width="60%">


# 3.实际应用

### 3.1 用于判断单个基因树的确定性，从而可以进行基因树的筛选
通过以1000次重复的bootstrap树为输入，最终可以通过计算单个基因树的TC值和TCA值来给不同的基因排序
我们可以选择确定性更高的树来进行建树期望得到更准确的结果

### 3.2 应用于性状
除了使用系统发育树划分，我们还可以用任意的可以将物种进行二分的性状数值来进行确定性的计算

### 3.3 利用 TC 和 TCA 评估数据分析中不同操作的影响
在进行一系列自认为能优化系统发育树的操作之后，比较基因树前后的TC 和TCA值的变化，可以用来确认操作是否有效

### 3.4 用于判断次要拓扑的分布情况
即IC ICA值的最基础用法

# 4. 缺陷
文中认为该方法的缺陷为无法判断低IC ICA值的原因是由于信息量的不足还是存在真实的系统发育冲突现象。
另外当基因树取样不全时无法使用RAxML进行IC ICA值的计算，计划于之后进行升级补充。但是目前的尝试发现是能够计算出具体数值的，因此可以认为此缺陷已经完成了修复

另外根据本人自己的实际操作，可以发现：
当我们是用QS 和IC值来反映同一个分支节点的时候，有时会出现结论相反的情况，例如QC值高但是IC值缺很低。目前可以认为这是二分法与quartet方法之间的不同而导致的，理论上划分更加细致的QS方法的结论可能是更加接近真实情况的。
