---
title: "Whole-genome sequences of Malawi cichlids reveal multiple radiations interconnected by gene flow"
authors: Milan Malinsky, Hannes Svardal, Alexandra M. Tyers, Eric A. Miska, Martin J. Genner, George F. Turner, Richard Durbin
year: 2018
citekey: malinskyWholegenomeSequencesMalawi2018
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Whole-genome sequences of Malawi cichlids reveal multiple radiations interconnected by gene flow
</div>

**Authors:** Milan Malinsky, Hannes Svardal, Alexandra M. Tyers, Eric A. Miska, Martin J. Genner, George F. Turner, Richard Durbin  
**Year:** 2018  
**Zotero:** [Open in Zotero](zotero://select/items/@malinskyWholegenomeSequencesMalawi2018)

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
  The hundreds of cichlid fish species in Lake Malawi constitute the most extensive recent vertebrate adaptive radiation. Here we characterize its genomic diversity by sequencing 134 individuals covering 73 species across all major lineages. The average sequence divergence between species pairs is only 0.1–0.25%. These divergence values overlap diversity within species, with 82% of heterozygosity shared between species. Phylogenetic analyses suggest that diversification initially proceeded by serial branching from a generalist Astatotilapia-like ancestor. However, no single species tree adequately represents all species relationships, with evidence for substantial gene flow at multiple times. Common signatures of selection on visual and oxygen transport genes shared by distantly related deep-water species point to both adaptive introgression and independent selection. These findings enhance our understanding of genomic processes underlying rapid species diversification, and provide a platform for future genetic analysis of the Malawi radiation.
</div>

---

# 1. 文献研究类群与使用的数据集
- 马拉维湖的马拉维湖慈鲷(*Malawi cichlids*)
- 采集了主要谱系的73个物种，134个个体。涉及到7个主要的生态类群  
- 使用全基因组测序数据
<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260415100904143.png width="60%">



# 2. 文献使用的主要方法
方法过多，不一一列举，主要集中在比较基因组学+系统发育学+分子生物学，又需要了解的部分，在下列具体内容中会指出。

# 3. 文献的主要内容与理论
##### **核心问题**：马拉维湖慈鲷经历了多次快速辐射演化，导致常规的线粒体基因和少量核基因的研究往往得出了相互矛盾的结论。因此亟需从全基因组层面来解答：这些物种到底是如何关联的？驱动它们快速适应极端环境（如深水区）的基因组机制是什么？

## 3.1 基因组差异性比较
具体而言，为了评估**物种内**的核苷酸多样性(nucleotide diversity (π))，测量了每个个体的**杂合位点频率**。同时也计算了**物种间平均两两序列差异(dxy)**。并进行了比较，得到如下结论：
- 马拉维湖慈鲷体内的基因变异水平即π，在动物界处于地位，说明其遗传背景较为纯净，彼此间应该较为相似，但其形态却千差万别  
- π值大小与种群的物种数量无相关性，说明**其新物种的可能更多受到生态位竞争或性选择的影响，而不是仅仅取决于遗传多样性的多寡**
- π值大小分布部分时候高于dxy，说明在部分极端情况下，**单个二倍体内部的序列差异，可能要高于两个不同物种之间的差异**。这极大地证明了这些物种的分化时间非常短
- 除了差异与多样性的比值较低外，大部分遗传变异在物种之间是共享的。平均而言，一个个体内部 82% 的杂合位点，其两个等位基因都能在其他物种中观察到，这与预期及此前观察到的高水平“不完全谱系分选”（ILS）相一致。
<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260415101722852.png width="60%">

## 3.2 较低的每代基因突变率
通过”家系测序“发现慈稠父代和子代之间的基因突变率仅为人类的1/3到1/4，但考虑到其1-2年一代的迭代频率，年平均突变率可能较高。综合而言，其世代突变速率并没有此前研究预计的哪些呈现一个极高的水平。说明**其物种的多样化可能并不是由于高世代突变率带来的**

## 3.3 基因组数据支持生态形态分类
根据其基因组变异数据（SNP），能够将马拉维湖的鱼类分为7类，与形态证据相符合。除了少数混和形态的鱼类以外。
- **Utaka 组**，其中一些物种与深海底栖组聚类更近，另一些则与浅海底栖组更近，**显示出基因流的可能**
- **孔雀鲷属（_Aulonocara_）** 的两个物种——_A. stuartgranti_ 和 _A. steveni_，它们位于浅海和深海底栖组之间。尽管它们像许多深海底栖物种（包括其他 _Aulonocara_ 属鱼类）一样拥有扩大的侧线感觉器官，但它们通常生活在较浅的水域。
<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260416153614399.png width="60%">  
## 3.4 等位基因共享与树状关系的不一致
根据上述结论，我们可以发现某些物种的位置十分摸棱两可，说明在基因组水平上可能处于明确定义的各个组别之间，说明不存在单一的进化树可以将这些物种联系起来
同时使用D统计量中的$D_{min}$来进行验证，在不预先设置任何亲缘关系的前提下，计算了$D_{min}$值，发现**62% 的三物种组合**（121,485 组中的 75,616 组）具有显著的正 $D_{min}$ 得分，揭示了多个层面的网状进化

## 3.5 系统发育树的构建
利用getWGSeq软件获得可用的全基因组多序列比对文件。设置了2543个无重叠的滑动窗口，每个窗口包含8000个SNP数据。
通过RAxML计算每个窗口的系统发育树，得到了2542个不同的系统发育树。随后使用诸多对ILS具有鲁棒性的方法合并物种树。例如：
**Bayesian SNAPP, SVDquartets, ASTRAL**等。并移除了此前基因组数据鉴定的“中间类群”，最终得到了较为一致性的系统发育物种树。

同时也构建了线粒体的物种树，最终得到了不一致的线粒体物种树，作者认为：
- 线粒体在ILS存在的情况下本就不足以反映正确的系统发育共识关系
- 线粒体本身也会收到自然选择的干扰
以此突出了使用全基因组数据的必要性
<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260416163204502.png width="60%">  

## 3.6 物种间渐渗信号的检测
将用于构建邻接树的成对遗传距离与样本沿树枝距离进行了对比，并计算了**残差（Residuals）**。发现了诸多差异，说明系统发育树无法解释物种间可能发生的杂交渐渗：
- **_P. cf. longimanus_**：它与深海底栖支系以及一部分浅海底栖支系（主要是 _Lethrinops_ 属）的遗传距离比进化树所暗示的更近。
- **_O. tetrastigma_（来自伊兰巴湖）**：它与广适性的 _A. calliptera_（特别是来自金吉里湖的样本，距离仅 3.2 公里）的距离比预期近得多。

此外发现了原本树中亲缘关系较远的物种之间共享了**长单倍型** ，此即近期发生了杂交渐渗的迹象

除此之外，还计算了根据ABBA统计量的计算原理设计了Fb统计量，来查看具体是哪些物种之间发生了显著的杂交渐渗。  
同时考虑，如果物种内存在明显的种群结构，加上物种形成事件快速连续发生，也可能导致fb值的显著升高。
然而，在被多次物种形成事件隔开的“远缘”支系之间存在多余的等位基因共享时，祖先种群结构必须穿过这些物种形成事件进行隔离，且不影响其姐妹支系，这种情况通常是不具有说服力的。因此，认为有强有力的证据表明存在多次跨物种的基因流事件。
<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260416164026855.png width="60%"> 

## 3.7 辐射演化的起源问题
_**A. calliptera**_ 在此前一直被认为是祖先态物种，但是全基因组显示其明确的作为岩栖组（mbuna）的姐妹群，仅有5.99%的滑动窗口树支持其作为祖先。
且在尝试了采样其它地区的 **_A. calliptera_** 后，其仍能聚成一个单系群，并在去除岩栖组（mbuna）后，仍保持系统发育位置不变，说明马拉维湖的慈稠嵌套位置并非由于后期的杂交所致。
因此推测其祖先很可能是一种在生态位和表型上与 A. calliptera 及其他 Astatotilapia 属物种相似的河流广适性鱼类。因为观察到其基因组上存在一定差异，但是在形态上存在相似性，**_A. calliptera_** 嵌套在其他亲缘关系较远但生态习性相似的 **_Astatotilapia_** 物种的形态空间内
<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260416173542868.png width="60%"> 

<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260416173558511.png width="60%"> 

因此作者提出了一个渐进式的三波爆发模型：
- 远洋组（pelagic）辐射最早被播种
- 接着是底栖组 + Utaka 组
- 最后是岩栖组 mbuna
将我们的世代突变率应用于观测到的基因组歧异度，我们得出这些支系之间的平均分化时间估计在 **46 万年前（ka）**至 **39 万年前**之间（假设每代三年）。这些估值点均落在马拉维湖古生态记录推测出的倒数第二次长期深水湖阶段内，而置信区间的上限则覆盖了约 80 万年前的第三次深水湖阶段。然而现存的所有 _A. calliptera_ 的共同祖先都非常年轻。
因此作者认为，产生这种现象的原因是，周边流域的_A. calliptera_在河流干涸期的时候大量灭绝，而马拉维湖因为蓄水量大所以大量的A. calliptera存活了下来。而后在涝期溢出到其它河流中。所以追溯到的共同祖先其实是最近的一次由马拉维湖流入到其它河流中的_A. calliptera_。

## 3.8 蛋白质编码基因的选择压力（寻找物种间存在差异的分子学原因）
定义了：
- **$p_S$ (同义突变频率)**：这种突变改变了 DNA 字母，但**没改变**氨基酸（蛋白质的成分）。它像背景噪音，不影响功能。
- **$p_N$ (非同义突变频率)**：这种突变**改变了**蛋白质的成分。如果这个改变对生存有利，它就会被自然选择保留下来
- 定义了受选择信号$Δ_{N-S} = $p_N$ - $p_S$
发现异常在于编码区的同义突变 $p_S$ 居然比非编码区高 13%。说明调控表型的aa序列比背景的相似度还要高。意味这可能是由于近期的杂交渐渗铺平了编码区间的差异

作者比较了慈鲷（_M. zebra_）和它的远亲斑马鱼（Zebrafish）的基因拷贝数，发现：
在两者的多拷贝基因家族中，受选择信号最高，这意味着导致慈稠多样性的原因可能并不是其特有的基因，而是来源于古老的多拷贝家族。并且通过GO富集鉴定了这些基因家族的功能，主要集中在：
- 血红蛋白功能和氧气运输
- 光传导与视觉感知
- 免疫系统

## 3.9 深度适应的共享机制
通过全基因组证据发现，亲缘关系较远的深水区鱼类之间，在深水适应性的基因上具有相似性。且发现$Δ_{N-S}$ 与Fdm高度相关，说明了这种深海适应性相关的基因似乎来源于渐渗。
#### 不同鱼类为何演化出了及其相似的性状且基因相似度也很高
文中给出了三种可能的情况以及对应的案例：
##### 1） 祖传变异的平行选择( Parallel Selection on Ancestral Variation )
祖先物种中留下了丰富的基因库，在后续分化出的不同物种中均保留了这一适应性的基因

##### 2） 适应性渗入(Adaptive Introgression)
通过杂交渐渗把实现不同物种间适应性的基因交流，例如：绿光敏感视蛋白基因（RH2 cluster），其在基因组中所处的位置检测出了明显的$f_{dM}$ 峰值    

<img src=https://github.com/XiongTor/Rosa_family/blob/main/paper/2026.04.27/file-20260416203007885.png width="60%"> 


##### 3）对新突变的独立选择 (Independent Selection on De Novo Mutations)---趋同演化
**外周蛋白-2（Peripherin-2）**和**抑制蛋白-C（Arrestin-C）**
$f_{dM}$ 的峰值非常窄，精准地在基因边界处消失。如果是杂交，周边的区域也会被带过来；既然只有基因内部有信号，说明不是杂交。


## 4. 思考
文章内容十分之全面
从比较基因组学入手，发现了马拉维湖鱼类基因组普遍较为纯净，变异较少，迭代变异速率较小，说明其多样并非由于高变异度导致，而是由于其它原因，结合后文应该是由于特性位置的部分基因的变化导致了其多样化。同时也暗示了存在杂交渐渗的可能，相互漂白了基因组。
其次从系统发育学与形态学的角度出发，结合比较基因组结果和系统发育树结果，推断出马拉维湖的鱼类存在复杂的网状进化事件，并提出了新的模型来解释马拉维湖的鱼类起源于演化历程。
最终落实到具体的基因与蛋白，明确了其调控通路。通过基因组进一步分析了深水区鱼类为何获得相同的性状。可能的原因各不相同。

总体而言，其研究范围不在于特定的类群，而是聚焦于马拉维湖这一区域，从生态位进化种群的划分，结合基因组与系统发育树深入探讨了其复杂的演化历史以及“一波三折”的辐射演化过程。最终落实到性状与基因，蛋白，详细验证落实了其分子层面的调控机制以及性状的获得机制。

