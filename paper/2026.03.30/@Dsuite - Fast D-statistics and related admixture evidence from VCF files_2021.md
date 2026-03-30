---
title: "Dsuite - Fast D-statistics and related admixture evidence from VCF files"
authors: Milan Malinsky, Michael Matschiner, Hannes Svardal
year: 2021
citekey: malinskyDsuiteFastDstatistics2021
tags: [paper, literature]
---

<div style="font-size: 28px; color: #C97C7C; margin-top: 0px;">
  Dsuite - Fast D-statistics and related admixture evidence from VCF files
</div>

**Authors:** Milan Malinsky, Michael Matschiner, Hannes Svardal  
**Year:** 2021  
**Zotero:** [Open in Zotero](zotero://select/items/@malinskyDsuiteFastDstatistics2021)

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
  Patterson's D, also known as the ABBA-BABA statistic, and related statistics such as the f4-ratio, are commonly used to assess evidence of gene flow between populations or closely related species. Currently available implementations often require custom file formats, implement only small subsets of the available statistics, and are impractical to evaluate all gene flow hypotheses across data sets with many populations or species due to computational inefficiencies. Here, we present a new software package Dsuite, an efficient implementation allowing genome scale calculations of the D and f4-ratio statistics across all combinations of tens or hundreds of populations or species directly from a variant call format (VCF) file. Our program also implements statistics suited for application to genomic windows, providing evidence of whether introgression is confined to specific loci, and it can also aid in interpretation of a system of f4-ratio results with the use of the “f-branch” method. Dsuite is available at https://github.com/millanek/Dsuite, is straightforward to use, substantially more computationally efficient than comparable programs, and provides a convenient suite of tools and statistics, including some not previously available in any software package. Thus, Dsuite facilitates the assessment of evidence for gene flow, especially across larger genomic data sets.
</div>

---

# 1. 基本原理
D检验原先是为群体遗传学服务，用于检测不同群体之间的基因交流现象。但同时也可以服务于研究近缘物种间的杂交与渗入，只需要满足以下群体遗传学的假设：
1. 物种间由于共同祖先和不完全谱系分选（ILS）而共享大量遗传变异 ；
2. 同一位点的回复突变和重复突变可以忽略不计 ；
3. 各物种间的替代速率是均一的

计算方法上，基于VCF文件，对所有可能的quartet组合计算其D值，即ABBA值
![](../imag/@Dsuite%20-%20Fast%20D-statistics%20and%20related%20admixture%20evidence%20from%20VCF%20files_2021/file-20260329132707149.png)

具体的计算公式如下：
其中p代表某一位点的等位基因频率
例如：群体 P1 有 5 个二倍体个体。在位点 $i$，有 3 个个体是 `0/1`，2 个个体是 `1/1`。
'1' 的总数 = $3 \times 1 + 2 \times 2 = 7$。
等位基因总数 = $5 \times 2 = 10$。
则 $p_{i1} = 7/10 = 0.7$。
![](../imag/@Dsuite%20-%20Fast%20D-statistics%20and%20related%20admixture%20evidence%20from%20VCF%20files_2021/file-20260329132811758.png)
计算D值
$$D = nABBA − nBABA / nABBA + nBABA$$


# 2. 使用事项
### 2.1 基础全局统计量，用于检测物种间的杂交渐渗现象
- D值计算
- f4值计算，排除了`ABBA  BABA`总数过少导致的D值偏大的系统误差
- f_branch值计算，更加准确的指出，渐渗是来源于近期的物种与物种之间，还是来源于早期的祖先杂交

### 2.2 滑动窗口统计量，用于检测特定染色体区域的渐渗基因------仅适用于物理上连续的基因片段
- fdm值计算，通过滑动窗口，计算fdm值定位渐渗基因可能发生的位置
- Df值计算，引入了 dxy，计算群体间的遗传距离，如果是近期杂交渐渗，dxy应该偏小，否则偏大
> 由于其本质上是按照等量的SNP来进行滑动窗口的计算，所以针对于物理上离散的，在染色体上并不连续的基因来说，该方法并不适用。因为当这些离散基因拼接成超矩阵时，实际上在染色体上并不一定是按照超矩阵的顺序排列，这意味着我们的离散基因其实可以随机排列，这样会导致每一个窗口的SNP并不是固定的，例如，第一次按照顺序1组成超矩阵，第二次按照顺序2组成超矩阵，会得到两个完全不同的结果

# 3.与其它软件的对比
- DMIXTOOLS
- HyDe
- comp-D
- PopGenome

Dsuite在运行速度和内存占用上均占优

# 4. 结果解读
主要是针对f_branch的结果解读：
- 图中的红色指针，是使用模拟数据得到的，一般情况下不会有该结果
- 图中是可以看出渐渗方向的，认为X轴的物种为P3，Y轴的物种为P2
![](../imag/@Dsuite%20-%20Fast%20D-statistics%20and%20related%20admixture%20evidence%20from%20VCF%20files_2021/file-20260329140616779.png)


# 5. 思考
- 否定了之前通过在离散的直系同源基因超矩阵上使用滑动窗口定位渐渗基因的想法，该做法太容易带来极大的不稳定性
- 明确了如何解读f_branch的结果，可以帮助确定一些高可能性的渐渗事件
- 再次明确了绝对遗传距离dxy可以帮助确定是近期渐渗还是其它可能




