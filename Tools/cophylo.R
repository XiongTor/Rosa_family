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

tree1 <- read.tree("tree/rosa_chloroplast_supermatrix.fasta.rt.tre")
tree1 <- ladderize(tree1, right = T)

tree2 <- read.tree("tree/rosa_ags353_treeshrink_sp_rt.tre")
tree2 <- ladderize(tree2, right = T)

tree3 <- read.tree("tree/rosa_orthofinder_MO_treeshrink_sp_rt_rename.tre")
tree3 <- ladderize(tree3, right = T)

tree4 <- read.tree("tree/rosa_ags353_concatenation_sp_partition.txt.rt.tre")
tree4 <- ladderize(tree4, right = T)

tree5 <- read.tree("tree/rosa_orthofinder_connect.txt.rt.tre")
tree5 <- ladderize(tree5, right = T)


# read the subfamily information
subf <- read.csv("Rosaceae_genus_accepted_lastest.csv")

# group
tips <- as.data.frame(tree1$tip.label)
colnames(tips) <- "species"
tips$genus <- sapply(strsplit(tips$species,"_"),`[`,1)

groups <- merge(tips,subf,by.x = "genus",by.y = "Genera",all.x = T) %>% 
  select(Subfamilies,Tribes,genus,species)

groups$Subfamilies[is.na(groups$Subfamilies)] <- "outgroups"

groups <- groups %>% 
  mutate(species=factor(species,levels = sort(tips$species))) %>% 
  arrange(species)

group_colors <- setNames(c("#4766b0","#CB793A","#006400","#A0A5A2"), unique(groups$Subfamilies))

link_colors <- make.transparent(group_colors[groups$Subfamilies], 0.4)

# count RF and path distance
rf_value_12 <- RF.dist(tree1, tree2, check.labels = TRUE,normalize = TRUE)
path_value_12 <- path.dist(tree1, tree2, check.labels = TRUE, use.weight = FALSE)

rf_value_13 <- RF.dist(tree1, tree3, check.labels = TRUE,normalize = TRUE)
path_value_13 <- path.dist(tree1, tree3, check.labels = TRUE, use.weight = FALSE)

rf_value_14 <- RF.dist(tree1, tree4, check.labels = TRUE,normalize = TRUE)
path_value_14 <- path.dist(tree1, tree4, check.labels = TRUE, use.weight = FALSE)

rf_value_15 <- RF.dist(tree1, tree5, check.labels = TRUE,normalize = TRUE)
path_value_15 <- path.dist(tree1, tree5, check.labels = TRUE, use.weight = FALSE)

rf_value_23 <- RF.dist(tree2, tree3, check.labels = TRUE,normalize = TRUE)
path_value_23 <- path.dist(tree2, tree3, check.labels = TRUE, use.weight = FALSE)

rf_value_45 <- RF.dist(tree4, tree5, check.labels = TRUE,normalize = TRUE)
path_value_45 <- path.dist(tree4, tree5, check.labels = TRUE, use.weight = FALSE)

rf_value_24 <- RF.dist(tree2, tree4, check.labels = TRUE,normalize = TRUE)
path_value_24 <- path.dist(tree2, tree4, check.labels = TRUE, use.weight = FALSE)

rf_value_35 <- RF.dist(tree3, tree5, check.labels = TRUE,normalize = TRUE)
path_value_35 <- path.dist(tree3, tree5, check.labels = TRUE, use.weight = FALSE)


# cophylo

jj_12 <- as.matrix(cbind(sort(tree1$tip.label),
                         sort(tree2$tip.label)))

jj_13 <- as.matrix(cbind(sort(tree1$tip.label),
                         sort(tree3$tip.label)))

jj_14 <- as.matrix(cbind(sort(tree1$tip.label),
                         sort(tree4$tip.label)))

jj_15 <- as.matrix(cbind(sort(tree1$tip.label),
                         sort(tree5$tip.label)))

jj_23 <- as.matrix(cbind(sort(tree2$tip.label),
                         sort(tree3$tip.label)))

jj_45 <- as.matrix(cbind(sort(tree4$tip.label),
                         sort(tree5$tip.label)))

jj_24 <- as.matrix(cbind(sort(tree2$tip.label),
                         sort(tree4$tip.label)))

jj_35 <- as.matrix(cbind(sort(tree3$tip.label),
                         sort(tree5$tip.label)))	

obj_12<- cophylo(tree1,tree2,assoc=jj_12,rotate=T)

obj_13<- cophylo(tree1,tree3,assoc=jj_13,rotate=T)

obj_14<- cophylo(tree1,tree4,assoc=jj_14,rotate=T)

obj_15<- cophylo(tree1,tree5,assoc=jj_15,rotate=T)

obj_23<- cophylo(tree2,tree3,assoc=jj_23,rotate=T)

obj_45<- cophylo(tree4,tree5,assoc=jj_45,rotate=T)

obj_24<- cophylo(tree2,tree4,assoc=jj_24,rotate=T)

obj_35<- cophylo(tree3,tree5,assoc=jj_35,rotate=T)


pdf("cytonuclear_final.pdf",width = 25,height = 20)
par(mfrow=c(2,4))  # 增加左侧边距为 8 行
par(oma = c(2, 2, 1.5, 2))
plot(obj_12,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("Chloroplast", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("AGS353_Coal", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_12), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_12), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_13,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("Chloroplast", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("Orthofinder_Coal", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_13), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_13), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_14,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("Chloroplast", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("AGS353_Concat", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_14), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_14), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_15,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("Chloroplast", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("Orthofinder_Concat", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_15), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_15), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_23,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("AGS353_Coal", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("Orthofinder_Coal", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_23), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_23), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_45,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("AGS353_Concat", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("Orthofinder_Concat", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_45), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_45), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_24,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("AGS353_Coal", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("AGS353_Concat", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_24), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_24), cex=1, side=1, line=0.5,at=0,col="black")

plot(obj_35,
     link.type = "curved",
     link.lwd = 3,
     link.lty = "solid",
     link.col = link_colors,
     fsize = 0.8)
mtext("Orthofinder_Coal", cex=1.2, side=2, line=-1.1,at=0.6,col="black")
mtext("Orthofinder_Concat", cex=1.2, side=4, line=-1.4,at=0.6,col="black")
mtext(paste0("RF distance: ",rf_value_35), cex=1, side=1, line=-1.1,at=0,col="black")
mtext(paste0("Path distance: ",path_value_35), cex=1, side=1, line=0.5,at=0,col="black")


dev.off()
