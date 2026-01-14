#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-12-18
# Description: Plot boxplot

# ==== Main Script Start ====
ggplot(data1,aes(x=Sequencing.strategies,y=gene.nuber))+
  geom_boxplot(aes(fill= Sequencing.strategies),width=0.35,linewidth= 0.8,outlier.shape = NA)+
  scale_fill_manual(values = col1)+
  geom_point(position = position_jitter(width = 0.1,height = 0))+
  scale_y_continuous(expand = c(0, 0),limits =c(0,360))+
  theme_classic()+
  labs(x="Sequencing strategies",y="Angiosperms353 count",fill="Sequencing strategies")+
  theme(text=element_text(family="serif"),
        legend.position = "none",
        #legend.position = c(0.13,0.85),
        #legend.key.size = unit(0.6, "lines"),
        #legend.spacing.y = unit(0.3, "lines"),
        #legend.text = element_text(size = 10),
        #legend.title = element_text(size = 10,face = "bold"),
        #legend.key.size = unit(0.8, "lines"),
       # legend.title = element_text(size = 10,face = "bold"),
        axis.line = element_line(size = 0.5,colour = "black"),
        axis.text = element_text(size = 19,colour = "black"),
        axis.title.y = element_text(vjust = 7,size=19),
        axis.title.x = element_text(vjust = -1,size=19),