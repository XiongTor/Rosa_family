# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.05.04
# This script is used to record the MSCquartets analysis
# install.packages("MSCquartets")

library(ape)
library(phangorn)
library(MSCquartets)
library(ggplot2)
library(tidyr)
library(dplyr)

### Reading the tree file and couning the quartets
sptree=read.tree("12-MSCquartet/rosa_orthofinder_MO_treeshrink_sp_rt_oneoutg_final.tre")
sptree2 <- write.tree(sptree)
sptree_list <- list(sptree)
class(sptree_list) <- "multiPhylo"

gtrees=read.tree("12-MSCquartet/rosa_orthofinder_MO_genetrees_oneoutg.tre")

### 得到数据集中所有物种名字
tnames=taxonNames(sptree_list)
   
# 统计基因树中所有可能的quartet的三种拓扑出现次数
QT=quartetTableParallel(gtrees,tnames,numCores=10)
#dim(QT)
save(QT,file = "QT.RData")

# 将QT中未分辨的计数分配给三种拓扑或者直接丢弃
RQT=quartetTableResolved(QT)
save(RQT,file = "RQT.RData")
### 根据共祖模型（multispecies coalescent, MSC）判断拓扑的分布情况是否与理论值显著不符
# 对应的α值
pTable=quartetTreeTestInd(RQT,model="T1",speciestree=sptree2)
pTable_2=quartetTreeTestInd(RQT,model="T3")

# 同上，但是是检验拓扑分布是否符合MSC star模型（1/3，1/3，1/3）预期的分布
# 对应的β值
Qtest=quartetStarTestInd(pTable)
Qtest_2=quartetStarTestInd(pTable_2)

#保存数据与再次读取数据
save(QT, RQT, pTable,pTable_2, Qtest,Qtest_2, file = "quartet_analysis_results.RData")


# 开始画图 ---------------------------------------------------------------------------------
# orthofinder
load("12-MSCquartet/orthofinder/quartet_analysis_results.RData")

# ags353
load("12-MSCquartet/ags353/quartet_analysis_results_3outg.RData")

#更改原始函数，方便画图调整：-----------------
simplexPrepare_scaled <- function(model = "T3", maintitle = NULL, titletext = NULL) {
  if (!(model %in% c("T1", "T3", "cut"))) 
    stop("Invalid model name; use 'T1','T3', or 'cut'.")
  
  lineWidth = 2
  oldpar = par(mar = c(0, 0, 4, 0) + 0.1) #调整图幅
  on.exit(par(oldpar))
  
  top = c(1, 0, 0)
  left = c(0, 1, 0)
  right = c(0, 0, 1)
  mid = rep(1, 3)
  
  # 在这里修改scale和shift的值来调整三角形
  scale = 0.5        # 改这个值调整大小，越小越大，此值为标尺
  shift_x = 0        # 改这个值左右移动
  shift_y = 0        # 改这个值上下移动
  
  lim <- c(-1, 1) * scale + shift_x
  yylim <- c(-0.1, 0.6) * scale + shift_y #调整图幅
  
  plot(0, 0, xlim = lim, ylim = yylim, type = "n", asp = 1, 
       axes = FALSE, xlab = "", ylab = "", main = maintitle,
       cex.main = 2) #调整主标题的大小
  
  mtext(eval(bquote(.(titletext))), side = 3, line = -1, cex = 1.4) #调整副标题大小
  
  MSCquartets:::simplexSegment(top, left, lty = "dotted", lwd = lineWidth)
  MSCquartets:::simplexSegment(top, right, lty = "dotted", lwd = lineWidth)
  MSCquartets:::simplexSegment(left, right, lty = "dotted", lwd = lineWidth)
  MSCquartets:::simplexSegment(top, mid, lwd = lineWidth)
  
  if (model == "T3") {
    MSCquartets:::simplexSegment(left, mid, lwd = lineWidth)
    MSCquartets:::simplexSegment(right, mid, lwd = lineWidth)
  }
  
  if (model == "cut") {
    botmid = c(0, 0.5, 0.5)
    rightmid = c(0.5, 0, 0.5)
    leftmid = c(0.5, 0.5, 0)
    MSCquartets:::simplexSegment(mid, botmid, lwd = lineWidth)
    MSCquartets:::simplexSegment(left, rightmid, lwd = lineWidth)
    MSCquartets:::simplexSegment(right, leftmid, lwd = lineWidth)
  }
}

# 步骤2: 临时替换
assignInNamespace("simplexPrepare", simplexPrepare_scaled, ns = "MSCquartets")

# 步骤4: 恢复原函数（可选，如果不需要可以不执行）
# assignInNamespace("simplexPrepare", MSCquartets:::simplexPrepare, ns = "MSCquartets")


#封装另一个函数，用于加上各个颜色点的数量占比
plot_with_stats <- function(pTable, test, alpha = 0.01, beta = 0.05, cex = 0.5,
                            legend_cex = 1.5, legend_pt_cex = 2,show_legend = TRUE) {
  # 复制原函数
  quartetTestPlot_modified <- MSCquartets::quartetTestPlot
  
  # 获取函数体
  body_code <- deparse(body(quartetTestPlot_modified))
  if (show_legend) {
    # 显示图例：正常替换参数
    body_code <- gsub('cex = 0.8', 'cex = 1.5', body_code)
    body_code <- gsub('pt.cex = 1.2', 'pt.cex = 2', body_code)
    body_code <- gsub('lty = 0, lwd = 2\\)', 'lty = 0, lwd = 2, bty = "n")', body_code)
  } else {
    # 不显示图例：直接删除包含legend的那几行
    # 找到legend开始到结束的所有行
    legend_start <- grep('legend\\("topleft"', body_code)
    if (length(legend_start) > 0) {
      # 删除legend那行及其后续的参数行
      body_code <- body_code[-c(legend_start:(legend_start+1))]
    }
  }
  
  # 重新构建函数
  body(quartetTestPlot_modified) <- parse(text = paste(body_code, collapse = "\n"))
  
  # 使用修改后的函数绘图
  quartetTestPlot_modified(pTable, test, alpha = alpha, beta = beta, cex = cex)
  
  # 计算统计
  p_test <- if(test == "T3") pTable[, "p_T3"] else pTable[, "p_T1"]
  p_star <- pTable[, "p_star"]
  p_test <- p_test[!is.na(p_test)]
  p_star <- p_star[!is.na(p_star)]
  n_total <- length(p_test)
  
  n_reds <- sum((p_test < alpha) & (p_star <= beta))
  n_blues <- sum((p_test >= alpha) & (p_star <= beta))
  n_yellows <- sum((p_test >= alpha) & (p_star > beta))
  n_oranges <- sum((p_test < alpha) & (p_star > beta))
  
  # 定义颜色
  redColor = "red"
  blueColor = "blue"
  yellowColor = "goldenrod3"
  orangeColor = "darkgreen"
  
  # 分别添加每种颜色的统计信息（两行显示）
  # 第一行：红色和蓝色
  text(x = -0.3, y = -0.35, 
       labels = sprintf("Red: %d (%.1f%%)", n_reds, n_reds/n_total*100),
       col = redColor, cex = 2, pos = 4)
  
  text(x = 0.1, y = -0.35, 
       labels = sprintf("Blue: %d (%.1f%%)", n_blues, n_blues/n_total*100),
       col = blueColor, cex = 2, pos = 4)
  
  # 第二行：黄色和绿色
  text(x = -0.3, y = -0.4, 
       labels = sprintf("Yellow: %d (%.1f%%)", n_yellows, n_yellows/n_total*100),
       col = yellowColor, cex = 2, pos = 4)
  
  text(x = 0.1, y = -0.4, 
       labels = sprintf("Green: %d (%.1f%%)", n_oranges, n_oranges/n_total*100),
       col = orangeColor, cex = 2, pos = 4)
}
# ---------------------------------------------------------------------------------------


# Qtest1, 抽取原来1/10的行
sample_size <- floor(nrow(Qtest) / 1000) # 计算 1/10 的行数，并向下取整
set.seed(123) # 为了结果的可重现性
sample_indices <- sample(1:nrow(Qtest),sample_size)
sampled_Qtest <- Qtest[sample_indices,]

# Qtest2, 抽取原来1/10的行
sample_size_2 <- floor(nrow(Qtest_2)/ 1000) # 计算 1/10 的行数，并向下取整
set.seed(123) # 为了结果的可重现性
sample_indices_2 <- sample(1:nrow(Qtest_2),sample_size_2)
sampled_Qtest_2 <- Qtest_2[sample_indices_2,]




### plot
alpha_list <- c(0.01,0.001,0.0001,0.00001,0.000001)
sample <- 0.001

for(i in alpha_list){
    pdf(paste0("./12-MSCquartet/ags353/MSCquartets_result_alpha",i,"_sample",sample,"_T3_sample_TEST.pdf"),width=16,height=18)
    par(mfrow=c(2,2))
    plot_with_stats(sampled_Qtest_2, "T3", alpha=i,beta = 0.05,cex = 1, 
                    legend_cex = 1.2, legend_pt_cex = 1.5) 
    plot_with_stats(sampled_Qtest_2, "T3", alpha=i,beta = 0.01, cex = 1, 
                    legend_cex = 1.2, legend_pt_cex = 1.5,show_legend = FALSE) 
    plot_with_stats(sampled_Qtest, "T1", alpha=i,beta = 0.05, cex = 1, 
                    legend_cex = 1.2, legend_pt_cex = 1.5,show_legend = FALSE) 
    plot_with_stats(sampled_Qtest, "T1", alpha=i,beta = 0.01, cex = 1, 
                    legend_cex = 1.2, legend_pt_cex = 1.5,show_legend = FALSE) 
    dev.off()
  }


## 计算并绘制两套不同的数据集，随着alpha值的变化，其红色点占比和黄色点占比的变化----------
alpha_list <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001)
beta_list <- c(0.01, 0.05)

# 初始化数据框
data_ags353 <- data.frame()

# 循环计算
for(alpha_val in alpha_list){
  for(beta_val in beta_list){
    # T1模型数据
    p_T1   <- Qtest[, "p_T1"]
    p_star <- Qtest[, "p_star"]
    
    # 去掉 NA
    p_T1   <- p_T1[!is.na(p_T1)]
    p_star <- p_star[!is.na(p_star)]
    
    # 总数
    n_total <- length(p_T1)
    
    # 各类统计（根据颜色分类）
    n_reds <- sum((p_T1 < alpha_val) & (p_star <= beta_val))        # reject tree & star
    n_blues <- sum((p_T1 >= alpha_val) & (p_star <= beta_val))      # fail to reject tree/reject star
    n_yellows <- sum((p_T1 >= alpha_val) & (p_star > beta_val))     # fail to reject tree & star
    n_oranges <- sum((p_T1 < alpha_val) & (p_star > beta_val))      # reject tree/fail to reject star
    
    # 比例
    red_ratio <- n_reds / n_total
    blue_ratio <- n_blues / n_total
    yellow_ratio <- n_yellows / n_total
    orange_ratio <- n_oranges / n_total
    
    # 保存T1数据
    data_ags353 <- rbind(data_ags353, data.frame(
      model = "T1",
      alpha = alpha_val,
      beta = beta_val,
      n_total = n_total,
      n_red = n_reds,
      n_blue = n_blues,
      n_yellow = n_yellows,
      n_orange = n_oranges,
      red_ratio = red_ratio,
      blue_ratio = blue_ratio,
      yellow_ratio = yellow_ratio,
      orange_ratio = orange_ratio
    ))
    
    # T3模型数据
    p_T3   <- Qtest_2[, "p_T3"]
    p_star <- Qtest_2[, "p_star"]
    
    # 去掉 NA
    p_T3   <- p_T3[!is.na(p_T3)]
    p_star <- p_star[!is.na(p_star)]
    
    # 总数
    n_total <- length(p_T3)
    
    # 各类统计
    n_reds <- sum((p_T3 < alpha_val) & (p_star <= beta_val))
    n_blues <- sum((p_T3 >= alpha_val) & (p_star <= beta_val))
    n_yellows <- sum((p_T3 >= alpha_val) & (p_star > beta_val))
    n_oranges <- sum((p_T3 < alpha_val) & (p_star > beta_val))
    
    # 比例
    red_ratio <- n_reds / n_total
    blue_ratio <- n_blues / n_total
    yellow_ratio <- n_yellows / n_total
    orange_ratio <- n_oranges / n_total
    
    # 保存T3数据
    data_ags353 <- rbind(data_ags353, data.frame(
      model = "T3",
      alpha = alpha_val,
      beta = beta_val,
      n_total = n_total,
      n_red = n_reds,
      n_blue = n_blues,
      n_yellow = n_yellows,
      n_orange = n_oranges,
      red_ratio = red_ratio,
      blue_ratio = blue_ratio,
      yellow_ratio = yellow_ratio,
      orange_ratio = orange_ratio
    ))
  }
}

data_ags353$type <- "ags353"
data_orthofinder$type <- "orthofinder"

# 合并数据
data_total <- rbind(data_ags353, data_orthofinder)

# 只保留T3模型的数据
data_total_T3 <- data_total %>%
  filter(model == "T3")

# 转换为长格式用于绘图
data_plot <- data_total_T3 %>%
  select(model, alpha, beta, red_ratio, yellow_ratio, type) %>%
  pivot_longer(cols = c(red_ratio, yellow_ratio), 
               names_to = "color_type", 
               values_to = "ratio")

# 计算y轴范围（统一两个图的范围）
y_max <- max(data_plot$ratio * 100)
y_min <- min(data_plot$ratio * 100)

# 绘制折线图
p_t3 <- ggplot(data_plot, aes(x = log10(alpha), y = ratio * 100, 
                      color = color_type, linetype = as.factor(beta))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ type, ncol = 2, 
             labeller = labeller(type = c("orthofinder" = "Orthofinder", 
                                          "ags353" = "AGS353"))) +
  scale_color_manual(values = c("red_ratio" = "red", "yellow_ratio" = "goldenrod3"),
                     labels = c("Red (reject both)", "Yellow (reject neither)")) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("β = 0.01", "β = 0.05")) +
  scale_y_continuous(limits = c(y_min, y_max)) +  # 统一y轴范围
  labs(x = "log10(alpha)", 
       y = "Percentage (%)",
       linetype = "Beta Value",
       title = "Red and Yellow Point Ratios vs Alpha (T3 Model)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 12, face = "bold"))

ggsave("./12-MSCquartet/Ratio_t3.pdf",p_t3,height = 4,width = 6)


### inferences the tree
#stree=QDC(gtrees)
#write.tree(stree,"QDS_tree.tre")

#stree=WQDC(gtrees)
#write.tree(stree,"WQDS_tree.tre")

### plot the network tree
#NANUQ(pTable, outfile = "NANUQdist",alpha = 0.000003, beta = 0.95, plot = FALSE)

