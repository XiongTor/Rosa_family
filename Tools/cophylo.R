#!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-12-18
# Description: Cophylo analysis

# ==== Main Script Start ====

### Compare topology
library("phytools")
library("dplyr")
library("ape")
library("phangorn")
library("ggplot2")

# In this example, we will compare the topology of three species trees: AGS353, orthofinder, and chloroplast.---2025.12.18
# read tree file 

ags353<- read.tree("cophylo/rosa_ags353_treeshrink_sp_rt.tre")
ags353 <- ladderize(ags353, right = T)

orthofinder_old <- read.tree("cophylo/rosa_orthofinder_sptree_rt.tre")
orthofinder_old <- ladderize(orthofinder_old, right = T)

orthofinder <- read.tree("cophylo/rosa_orthofinder_sp_rt_new.tre")
orthofinder <- ladderize(orthofinder, right = T)

chloroplast <- read.tree("cophylo/rosa_chloroplast_sp_rt.tre")
orthofinder <- ladderize(orthofinder, right = T)

# read the subfamily information
subf <- read.csv("cophylo/Rosaceae_genus_accepted_lastest.csv")

# group
tips <- as.data.frame(ags353$tip.label)
colnames(tips) <- "species"
tips$genus <- sapply(strsplit(tips$species,"_"),`[`,1)

groups <- merge(tips,subf,by.x = "genus",by.y = "Genera",all.x = T) %>% 
  select(genus,species,Subfamilies)

groups$Subfamilies[is.na(groups$Subfamilies)] <- "outgroups"

groups <- groups %>% 
  mutate(species=factor(species,levels = sort(tips$species))) %>% 
  arrange(species)

group_colors <- setNames(c("#4766b0","#CB793A","#006400","#A0A5A2"), unique(groups$Subfamilies))

link_colors <- make.transparent(group_colors[groups$Subfamilies], 0.4)

# count RF and path distance
rf_value_a_c <- RF.dist(ags353, chloroplast, check.labels = TRUE,normalize = TRUE)
path_value_a_c <- path.dist(ags353, chloroplast, check.labels = TRUE, use.weight = FALSE)

rf_value_o_c <- RF.dist(orthofinder, chloroplast, check.labels = TRUE,normalize = TRUE)
path_value_o_c <- path.dist(orthofinder, chloroplast, check.labels = TRUE, use.weight = FALSE)

rf_value_a_o <- RF.dist(ags353, orthofinder, check.labels = TRUE,normalize = TRUE)
path_value_a_o <- path.dist(ags353, orthofinder, check.labels = TRUE, use.weight = FALSE)


rf_value_o_o <- RF.dist(orthofinder_old, orthofinder, check.labels = TRUE,normalize = TRUE)
path_value_o_o <- path.dist(orthofinder_old, orthofinder, check.labels = TRUE, use.weight = FALSE)
# cophylo

jj_a_c <- as.matrix(cbind(sort(ags353$tip.label),
                      sort(chloroplast$tip.label)))

jj_o_c <- as.matrix(cbind(sort(orthofinder$tip.label),
                          sort(chloroplast$tip.label)))

jj_a_o <- as.matrix(cbind(sort(ags353$tip.label),
                          sort(orthofinder$tip.label)))

jj_o_o <- as.matrix(cbind(sort(orthofinder_old$tip.label),
                          sort(orthofinder$tip.label)))

obj_a_c<- cophylo(ags353,chloroplast,assoc=jj_a_c,rotate=T)

obj_o_c<- cophylo(orthofinder,chloroplast,assoc=jj_a_o,rotate=T)

obj_a_o<- cophylo(ags353,orthofinder,assoc=jj_a_o,rotate=T)

obj_o_o<- cophylo(orthofinder_old,orthofinder,assoc=jj_a_o,rotate=T)



pdf("cytonuclear.pdf",width = 15,height = 10)
par(mfrow=c(1,3))  # 增加左侧边距为 8 行
par(oma = c(2, 2, 1.5, 2))
plot(obj_a_c,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("AGS353", cex=1.2, side=2, line=-0.6,at=0.6,col="black")
mtext("Chloroplast", cex=1.2, side=4, line=-0.8,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_a_c), cex=1, side=1, line=-1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_a_c), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_o_c,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("orthofinder", cex=1.2, side=2, line=-1.6,at=0.6,col="black")
mtext("Chloroplast", cex=1.2, side=4, line=-0.8,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_o_c), cex=1, side=1, line=-1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_o_c), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_a_o,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("Ags353", cex=1.2, side=2, line=-1.3,at=0.6,col="black")
mtext("orthofinder", cex=1.2, side=4, line=-0.6,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_a_o), cex=1, side=1, line=-1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_a_o), cex=1, side=1, line=0.5,at=0,col="black")

dev.off()