#### Load packages ----
library(Seurat)
library(tidyverse)
library(ggpubr)

##
public <- readRDS("./2_anno.rds")
head(public@meta.data)
table(public$cell_type)

##
public_End <- subset(public, idents = "Endothelial")
matrix <- as.data.frame(t(as.data.frame(public_End@assays[["RNA"]]@data)))
expression <- matrix[c("PCAT19","MALAT1")]
meta <- public_End@meta.data[,c(1,21)]
data <- cbind(meta,expression)

data1 <- data[,c(2,3)]
data1$gene <- "PCAT19"
colnames(data1) <- c("group","expression","gene")
data11=data1[which(rowSums(data1==0)==0),]

data2 <- data[,c(2,4)]
data2$gene <- "MALAT1"
colnames(data2) <- c("group","expression","gene")
data22=data2[which(rowSums(data2==0)==0),]

data3 <- rbind(data11,data22)

data3$group[data3$group == "contol"] <- "ctrl"
data3$group[data3$group == "case"] <- "atherosclerosis"
p <- ggplot(data=data3, aes(x = group, y = expression,fill = group)) + 
  geom_violin(trim = T, position = position_dodge(width = 1),scale = "width",lwd=0.3,alpha=0.8)+
  geom_boxplot(position = position_dodge(1),width=.2,lwd=0.3,alpha=0.5,aes(fill = group))+
  geom_point(size=0,alpha=0) +
  scale_fill_manual(values = c("#FE8D3C", "#56B4E9"))+
  theme_bw()+
  theme(axis.text.x = element_text(colour="black",family="",size=8), 
        axis.text.y = element_text(family="",size=8,colour="black"), 
        axis.title.y = element_text(family="",size = 8), 
        axis.ticks.x = element_line(colour="black",size=0.2),
        axis.line.x = element_line(colour="black",size=0.2),
        axis.line.y = element_line(colour="black",size=0.2),
        strip.text.x = element_text(size = 7),  
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
        legend.text=element_text(family="", colour="black", size=8), 
        legend.title=element_text(family="", colour="black", size=8),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Expression")+xlab("")+ #
  facet_wrap(.~gene,scales = "free",nrow = 2)+
  NoLegend()
p
ggsave(file.path(output.dir, "End_Vlnplot_MALAT1_PCAT19.pdf"),
       plot=p, height=4, width=1.8)

##
public_SMC <- subset(public, idents = "Smooth muscle cell")
matrix <- as.data.frame(t(as.data.frame(public_SMC@assays[["RNA"]]@data)))
expression <- matrix[c("SNHG9","MIR100HG","SNHG29","NEAT1")]
meta <- public_SMC@meta.data[,c(1,21)]
data <- cbind(meta,expression)

data1 <- data[,c(2,3)]
data1$gene <- "SNHG9"
colnames(data1) <- c("group","expression","gene")
data11=data1[which(rowSums(data1==0)==0),]

data2 <- data[,c(2,4)]
data2$gene <- "MIR100HG"
colnames(data2) <- c("group","expression","gene")

data3 <- rbind(data11,data2)
data3$group[data3$group == "control"] <- "ctrl"
data3$group[data3$group == "case"] <- "atherosclerosis"
data3$gene <- factor(data3$gene, levels=c("SNHG9","MIR100HG"))
data3$group <- factor(data3$group, levels=c("atherosclerosis","ctrl"))
boxplot <- ggpaired(data3, 
                    x = "group", y ="expression", color="group",
                    # xlab=x, ylab=y,
                    # legend.title=Macrophage, 
                    # font.label = list(size = 20, color = "black"),
                    palette = c("#e23620","#3498DB"),#"jco",
                    point.size = 0.2,width = 0.5,
                    line.color = "gray", line.size = 0.4,
                    short.panel.labs = FALSE)+
  # geom_jitter(width = 0.2)+
  labs(title = "",
       y = "Expression")+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        # aspect.ratio = 1.5,
        strip.text.x = element_text(size = 7),  
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=0.5),
        plot.title = element_text(face = "italic",size=6),
        axis.title = element_text(size=7), 
        axis.text.x = element_text(size=6.3),
        axis.text.y = element_text(size=6.3),
        axis.title.x = element_blank()) +
  theme_cat() +
  ggsignif::geom_signif(comparisons = list(c("atherosclerosis", "ctrl")),
                        test = "wilcox.test",
                        y_position = 2.8,
                        map_signif_level = F, 
                        textsize = 2.2,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) + 
  facet_wrap(.~gene,scales = "free",nrow = 2)+
  NoLegend() 
boxplot
ggsave(file.path(output.dir, "SMC_Boxplot_SNHG9-MIR100HG.pdf"),
       plot = boxplot, height = 4, width = 1.3)


##
public_Mac <- subset(public, idents = "Macrophage")
matrix <- as.data.frame(t(as.data.frame(public_Mac@assays[["RNA"]]@data)))
expression <- matrix[c("FTX","SNHG8","NEAT1")]
meta <- public_Mac@meta.data[,c(1,21)]
data <- cbind(meta,expression)

data1 <- data[,c(2,3)]
data1$gene <- "FTX"
colnames(data1) <- c("group","expression","gene")
data11=data1[which(rowSums(data1==0)==0),]

data2 <- data[,c(2,4)]
data2$gene <- "SNHG8"
colnames(data2) <- c("group","expression","gene")
data22=data2[which(rowSums(data2==0)==0),]

data3 <- rbind(data11,data22)
data3$group[data3$group == "contol"] <- "ctrl"
data3$group[data3$group == "case"] <- "atherosclerosis"
data3$gene <- factor(data3$gene,levels=c("SNHG8","FTX"))
p <- ggplot(data=data3, aes(x = group, y = expression,fill = group)) + 
  geom_violin(trim = T, position = position_dodge(width = 1),scale = "width",lwd=0.3,alpha=0.8)+
  geom_boxplot(position = position_dodge(1),width=.2,lwd=0.3,alpha=0.5)+
  geom_point(size=0,alpha=0) +
  scale_fill_manual(values = c('#FDC9C8',"#877FB2"))+
  theme_bw()+
  theme(axis.text.x = element_text(colour="black",family="",size=7), 
        axis.text.y = element_text(family="",size=7,colour="black"),
        axis.title.y = element_text(family="",size = 7), 
        axis.ticks.x = element_line(colour="black",size=0.2),
        axis.line.x = element_line(colour="black",size=0.2),
        axis.line.y = element_line(colour="black",size=0.2),
        strip.text.x = element_text(size = 7), 
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(family="", colour="black", size=6), 
        legend.title=element_text(family="", colour="black", size=6),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Expression")+xlab("")+ 
  stat_compare_means(aes(group = group),
                     size=2,
                     #label = "..p.signif..",
                     label = "p.format",
                     method = "wilcox.test",
                     hide.ns = T)+
  facet_wrap(.~gene,scales = "free",nrow = 1)+NoLegend()
p
ggsave(file.path(output.dir, "Mac_Vlnplot_FTX_SNHG8.pdf"),
       plot = p, height=1.8, width = 3)
