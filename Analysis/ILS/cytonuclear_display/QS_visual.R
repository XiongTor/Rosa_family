# ！/usr/bin/R
# Author:Xiongtao 

#list.files(path = ".",pattern = ".csv",full.names = T)

suppressWarnings(suppressMessages(library(treeio)))
suppressWarnings(suppressMessages(library(ggtree)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ape)))

#get the path
path <- commandArgs(TRUE)

qc <- read.tree(paste0(path,"/","QS_result/RESULT.labeled.tre.qc"))
qd <- read.tree(paste0(path,"/","QS_result/RESULT.labeled.tre.qd"))
qi <- read.tree(paste0(path,"/","QS_result/RESULT.labeled.tre.qi"))


# process node labels of above three labeled trees
# qc tree
tree <- qc
tree$node.label <- gsub("qc=","",tree$node.label)
#View(tree$node.label)
write.tree(tree,paste0(path,"/","QS_visual/tree_qc.tre"))
# qd tree
tree <- qd
tree$node.label <- gsub("qd=","",tree$node.label)
#View(tree$node.label)
write.tree(tree,paste0(path,"/","QS_visual/tree_qd.tre"))
# qi tree
tree <- qi
tree$node.label <- gsub("qi=","",tree$node.label)
#View(tree$node.label)
write.tree(tree,paste0(path,"/","QS_visual/tree_qi.tre"))


# read 3 modified tree files for QC/QD/QI
tree_qc <- read.newick(paste0(path,"/","QS_visual/tree_qc.tre"), node.label='support')
tree_qd <- read.newick(paste0(path,"/","QS_visual/tree_qd.tre"), node.label='support')
tree_qi <- read.newick(paste0(path,"/","QS_visual/tree_qi.tre"), node.label='support')


# add a customized label for internode or inter-branch, i.e., qc/qd/qI
score_raw = paste(tree_qc@data$support,"/",tree_qd@data$support,"/",tree_qi@data$support,sep="")
score = gsub("NA/NA/NA","",score_raw)
score = gsub("NA","-",score)
#View(score)


# set labeled QC tree as the main plot tree
tree <- tree_qc
tree@data$score <- score


#####################################################
# Partitioning quartet concordance. QC values were divided into four categories and this
# information was used to color circle points. 

# drop the internodes without QC vaule
#root <- tree@data$node[is.na(tree@data$support)]

pdf(file="treeQC_rect_circ.pdf", width = 60, height = 60) 

# (1)color circle points
ggtree(tree,layout = "re",branch.length = "none",size=1) +
  geom_tiplab(size=7) +
  xlim(0,6.5)+
  #geom_text(aes(label=node),size=2,hjust=-1)+
  geom_nodepoint(aes(subset=!isTip, fill=cut(support, c(1, 0.2, 0, -0.05, -1))),
                 shape=21, size=8) +
  theme(legend.position=c(0.2, 0.9),
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18))+
  scale_fill_manual(values=c("#2F4F4F", "#98F898", "#FFCC99","#FF6600"),
                    guide="legend", name="Quartet Concordance(QC)",
                    breaks =c("(0.2,1]","(0,0.2]","(-0.05,0]","(-1,-0.05]"),
                    labels=expression(QC>0.2, 0 < QC * " <= 0.2", -0.05 < QC * " <= 0", QC <= -0.05))


# (2)color branch
ggtree(tree, aes(color=cut(support, c(1, 0.2, 0, -0.05, -1))), layout="circular", size=1.3) +
  theme(legend.position=c(0.2, 0.9),
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18)) +
  scale_colour_manual(values=c("#2F4F4F", "#98F898", "#FFCC99","#FF6600"),
                      breaks=c("(0.2,1]","(0,0.2]","(-0.05,0]","(-1,-0.05]"),
                      na.translate=T, na.value="gray",
                      guide="legend", name="Quartet Concordance(QC)",
                      labels=expression(QC>0.2, 0 < QC * " <= 0.2", -0.05 < QC * " <= 0", QC <= -0.05))
# (3)color circle points and label each internode with QC/QD/QI

ggtree(tree, size=0.5,branch.length = "none") +
  geom_tiplab(size=7) +
  xlim(0,6.5)+
  geom_nodepoint(aes(subset=!isTip, fill=cut(support, c(1, 0.2, 0, -0.05, -1))),
                 shape=21, size=8) +
  theme(legend.position=c(0.2, 0.9),
        legend.key.size = unit(35, "pt"),
        legend.text = element_text(size=18),
        legend.title = element_text(size=18)) +
  scale_fill_manual(values=c("#2F4F4F", "#98F898", "#FFCC99","#FF6600"),
                    guide="legend", name="Quartet Concordance(QC)",
                    breaks=c("(0.2,1]","(0,0.2]","(-0.05,0]","(-1,-0.05]"),
                    labels=expression(QC>0.2, 0 < QC * " <= 0.2", -0.05 < QC * " <= 0", QC <= -0.05))+
  geom_text(aes(label=score, x=branch), size=7,vjust=-0.6)

dev.off()
