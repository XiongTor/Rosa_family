#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-10-16
# Description: Comparison of D/f4 values between and within subfamilies 

# ==== Main Script Start ====


#usage:Rscript plot_box.R f4.tsv p.tsv sub_trib.csv
library(ape)
library(dplyr)
library(ggtree)
library(tidyverse)

file <- commandArgs(trailingOnly =TRUE)


f4 <- read.table(file[1])
f4 <- f4[!(colnames(f4) %in% c("Outgroup","Morus_indica","Elaeagnus_angustifolia")),
         !(rownames(f4) %in% c("Outgroup","Morus_indica","Elaeagnus_angustifolia"))]

p <- read.table(file[2])
p <- p[!(colnames(p) %in% c("Outgroup","Morus_indica","Elaeagnus_angustifolia")),
       !(rownames(p) %in% c("Outgroup","Morus_indica","Elaeagnus_angustifolia"))]
sub_trib <- read.csv(file[3])

# Rosoideae
extract_subfamily_data <- function(subfamily1,subfamily2,f4,p,sub_trib){
  
  data <- f4 [colnames(f4) %in% sub_trib$Accepted_name[sub_trib$Subfamily == subfamily1],
              row.names(f4) %in% sub_trib$Accepted_name[sub_trib$Subfamily == subfamily2]]
  
  data_long <- data %>% 
    rownames_to_column("row") %>% 
    pivot_longer(-row, names_to = "col", values_to = "f4") %>%
    filter(row < col) %>% 
    mutate(pair=paste(row,col,sep="_"),type=paste0(subfamily1,"_",subfamily2)) %>% 
    select(pair,f4,type)
  
  data_p <- p[colnames(p) %in% sub_trib$Accepted_name[sub_trib$Subfamily == subfamily1],
              row.names(p) %in% sub_trib$Accepted_name[sub_trib$Subfamily == subfamily2]]
  
  data_long_p <- data_p %>% 
    rownames_to_column("row") %>% 
    pivot_longer(-row, names_to = "col", values_to = "p") %>%
    filter(row < col) %>% 
    mutate(pair=paste(row,col,sep="_"),type=paste0(subfamily1,"_",subfamily2)) %>% 
    select(pair,p,type)
  
  data_final <- merge(data_long,data_long_p,by=c("pair","type"),all.x=T) %>% dplyr::filter(p < 0.001)
  
  return(data_final)
}

Rosoideae <- extract_subfamily_data("Rosoideae","Rosoideae",f4,p,sub_trib)

Dryadoideae <- extract_subfamily_data("Dryadoideae","Dryadoideae",f4,p,sub_trib)

Amygdaloideae <- extract_subfamily_data("Amygdaloideae","Amygdaloideae",f4,p,sub_trib)

Dryadoideae_Rosoideae <- extract_subfamily_data("Dryadoideae","Rosoideae",f4,p,sub_trib)

Rosoideae_Amygdaloideae <- extract_subfamily_data("Rosoideae","Amygdaloideae",f4,p,sub_trib)

Dryadoideae_Amygdaloideae <- extract_subfamily_data("Dryadoideae","Amygdaloideae",f4,p,sub_trib)
#plot 
final <- 
  rbind(Rosoideae,
        Dryadoideae,
        Amygdaloideae,
        Dryadoideae_Rosoideae,
        Rosoideae_Amygdaloideae,
        Dryadoideae_Amygdaloideae)

final$type <- factor(final$type, levels = c("Rosoideae_Rosoideae", "Amygdaloideae_Amygdaloideae", "Dryadoideae_Dryadoideae", "Dryadoideae_Rosoideae","Dryadoideae_Amygdaloideae","Rosoideae_Amygdaloideae"))

final$type <- recode(final$type,
                     "Rosoideae_Rosoideae" = "RR",
                     "Amygdaloideae_Amygdaloideae" = "AA",
                     "Dryadoideae_Dryadoideae" = "DD",
                     "Dryadoideae_Rosoideae" = "DR",
                     "Dryadoideae_Amygdaloideae" = "DA",
                     "Rosoideae_Amygdaloideae" = "RA")

summary_final <- final %>% 
  group_by(type) %>%
  summarise(median_f4 = median(f4),
            mean_f4 = mean(f4),
            sd_f4 = sd(f4),
            n=n(),
            se_f4=sd_f4/sqrt(n))


col1=c("#3E5641","#A24936","#AE6C25","#D36135","#E3B505","#BBBE64")

p1 <- ggplot(summary_final, aes(x = type, y = mean_f4, fill = type)) +
  geom_col(color = "black", width = 0.5) +  # 柱状图
  geom_errorbar(aes(ymin = mean_f4 - se_f4, ymax = mean_f4 + se_f4),
                width = 0.2, size = 0.8) +  # 误差线
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +  # 添加y=0的参考线
  labs(x = "Type", y = "Mean f4 ± SE") +
  theme_bw(base_size = 14) +
  scale_fill_manual(values = col1) +
  theme(
    axis.text.x = element_text(hjust = 1, vjust = 1),
    legend.position = "none"
  )


p2 <- ggplot(data=final,aes(x=type, y=f4, fill=type))+
  geom_boxplot(outlier.size = 0.5)+
  geom_jitter(shape=16, position=position_jitter(0.2), size=0.5, alpha=0.5)+
  theme_classic()+
  scale_fill_manual(values = col1) +
  theme(axis.text.x = element_text(hjust = 1),
        legend.position = "none")

ggsave("mean_subfamily.pdf",p1,width= 6 , height= 6)
ggsave("boxplot_subfamily.pdf",p2,width= 6 , height= 6)

