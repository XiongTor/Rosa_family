library(ggplot2)
library(dplyr)

# 1. 预处理：按组 (Sample_Name) 寻找每条线的最大值
df <- read.csv("extracted_phylogenetic_data.csv") %>%
  mutate(across(c(Sim_Value, R2_Value), as.numeric)) %>%
  arrange(Sim_Value) %>%
  # 转换为因子确保等宽
  mutate(Sim_Value = factor(Sim_Value, levels = unique(Sim_Value))) %>%
  # 关键：按组标注每条线的最大值点
  group_by(Sample_Name) %>%
  mutate(is_group_max = (R2_Value == max(R2_Value, na.rm = TRUE))) %>%
  ungroup()

# 2. 绘图
ggplot(df, aes(x = Sim_Value, y = R2_Value, color = Sample_Name, group = Sample_Name)) +
  geom_line(linewidth = 0.8) + # 线条略透明，突出点
  geom_point(size = 1.5) +    # 普通点缩小变淡
  # 绘制每条线自己的最大值：使用实心菱形 (shape=18) 或 实心五角星 (shape=11且加厚)
  geom_point(data = filter(df, is_group_max), 
             shape = 18,   # 18是实心菱形，科研图表常用，比11号星星扎实
             size = 4) +   # 增大尺寸，颜色会自动跟随分组
  theme_bw() +
  labs(x = "Sim Level", y = expression(R^2), title = "Per-group Max R2 Highlighted") +
  theme(panel.grid.minor = element_blank(), legend.position = "right")
