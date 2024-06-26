##
library(Seurat)
library(tidyverse)
library(ggplot2)
output.dir <- paste0("~/results-V23-Figure/Fig1") 

##
seurat.obj <- readRDS("./anno/harmony_anno_20_0.8.rds")

##
a <- read.csv("./results-V23-Figure/Fig1/lncRNA/4-final_lncRNA.csv")
a <- a$lncRNA

counts <- as.data.frame(seurat.obj[["RNA"]]@counts)
b <- rownames(counts)

sample_list <- list(c1 = b, c2 = a)

library(VennDiagram)
library(readr)
library(eulerr)
p <- plot(
  euler(sample_list, shape = "ellipse"),
  quantities = T,
  labels = c("mRNA", "lncRNA"),
  edges = list(lty = 1,
               lex = 0.3,
               col = 'black'),
  fill = c('#95c2e7','#f7b496','#f8a6ac')) 
ggsave(file.path(output.dir, "Vein_mRNA-lncRNA.pdf"),
       plot = p,
       height = 2.5,
       width = 4.5)

##
markers <- c("C1QA","C1QB","C1QC",
             "CCL14","ACKR1","PLVAP",
             "APOD","PLA2G2A","FBLN1",
             "TPM2","PLN","MYH11",
             "IL32","TRAC","IL7R",
             "BGN","FN1","IGFBP2",
             "FABP4","ID4","MUSTN1",
             "CD79A","CD37","LTB",
             "IGLC2","IGHM","JCHAIN",
             "NKG7","GNLY","PRF1",
             "S100B","PLP1","GPM6B",
             "TPSB2","TPSAB1","AREG")
gene_cell_exp1 <- AverageExpression(seurat.obj,
                                    features = markers,
                                    group.by = 'cell_type',
                                    slot = 'data')
gene_cell_exp1 <- as.data.frame(gene_cell_exp1$RNA)

library(ComplexHeatmap)
df1 <- data.frame(colnames(gene_cell_exp1))
colnames(df1) <- 'class'
top_anno1 = HeatmapAnnotation(df = df1,
                              border = T,
                              show_annotation_name = F,
                              gp = gpar(col = 'black'),
                              col = list(class = c("Macrophage"="#00CD9B",
                                                   "Endothelial"="#C3A5F7FF",
                                                   "Fibroblast"="#098B8B",
                                                   "SMC"="#F8766D",
                                                   'T cell'="#FFB8A4FF",
                                                   "Fibromyocyte"="#D8C6F3FF",
                                                   "Pericyte"='#B53E2B',
                                                   "B cell"="#ACD45E",
                                                   "Plasma cell"='#8C549C',
                                                   "NK cell"="#9ECABE",
                                                   "Neuron"="#2F528F",
                                                   "Mast cell"="#E3AD68")))
#数据标准化缩放一下
library(circlize)
marker_exp1 <- t(scale(t(gene_cell_exp1),scale = T,center = T))
ht1 <- Heatmap(marker_exp1,
               row_title = "",
               column_title = "Cell Type Enrichment Score - mRNA",
               column_title_gp = gpar(fontsize = 8, fontface = "bold"),
               cluster_rows = F,cluster_columns = F,
               row_names_side = 'right',row_names_gp = gpar(fontface = 'italic',fontsize = 8,col="black"), #设置基因字体大小
               column_names_side = 'top',column_names_rot = 45,column_names_gp = gpar(fontsize = 8,col="black"),
               heatmap_legend_param = list(title = "Row Z Score", 
                                           title_gp = gpar(fontsize = 8,col="black"),
                                           labels_gp = gpar(fontsize = 8,col="black")),
               border = 'black',
               rect_gp = gpar(col = "white", lwd = 1),
               top_annotation = top_anno1,
               col = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
)
pdf(file.path(output.dir, "heatmap_mRNA.pdf"),width=5.2,height=5.5)
print(ht1)
dev.off()

##
lncRNA <- c("FILNC1","LINC00543", "HCG24","PCAT19", "FENDRR", "LINC02668",  
            "LINCADL","LINC02518","FAM3D-AS1","LINC01615","CARMN","ADAMTS9-AS1",          
            "LINC00402","PCED1B-AS1","LOH12CR2","LINC01902", "FLG-AS1", "IL12A-AS1",       
            "LINC00924","H19", "LINC01142","LINC01781","LINC00926","TRERNA1",
            "KCNQ1-AS1","LINC01055", "FRGCA","PRKCQ-AS1", "LINC01871","LINC02332",  
            "CHL1-AS2","FOXD3-AS1","LINC01608","LINC00520","HCG20","ARHGAP31-AS1")   
gene_cell_exp2 <- AverageExpression(seurat.obj,
                                    features = lncRNA,
                                    group.by = 'cell_type',
                                    slot = 'data')
gene_cell_exp2 <- as.data.frame(gene_cell_exp2$RNA)

library(ComplexHeatmap)
df2 <- data.frame(colnames(gene_cell_exp2))
colnames(df2) <- 'class'
top_anno2 = HeatmapAnnotation(df = df2,
                              border = T,
                              show_annotation_name = F,
                              gp = gpar(col = 'black'),
                              col = list(class = c("Macrophage"="#00CD9B",
                                                   "Endothelial"="#C3A5F7FF",
                                                   "Fibroblast"="#098B8B",
                                                   "SMC"="#F8766D",
                                                   'T cell'="#FFB8A4FF",
                                                   "Fibromyocyte"="#D8C6F3FF",
                                                   "Pericyte"='#B53E2B',
                                                   "B cell"="#ACD45E",
                                                   "Plasma cell"='#8C549C',
                                                   "NK cell"="#9ECABE",
                                                   "Neuron"="#2F528F",
                                                   "Mast cell"="#E3AD68")))
#数据标准化缩放一下
marker_exp2 <- t(scale(t(gene_cell_exp2),scale = T,center = T))
library(circlize)
ht2 <- Heatmap(marker_exp2,
               row_title = "",
               column_title = "Cell Type Enrichment Score - lncRNA",
               column_title_gp = gpar(fontsize = 8, fontface = "bold"),
               cluster_rows = F,cluster_columns = F,
               row_names_side = 'right',row_names_gp = gpar(fontface = 'italic',fontsize = 8,col="black"), #设置基因字体大小
               column_names_side = 'top',column_names_rot = 45,column_names_gp = gpar(fontsize = 8,col="black"),
               heatmap_legend_param = list(title = "Row Z Score", 
                                           title_gp = gpar(fontsize = 8,col="black"),
                                           labels_gp = gpar(fontsize = 8,col="black")),
               border = 'black',
               rect_gp = gpar(col = "white", lwd = 1),
               top_annotation = top_anno2,
               col = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033")))
ht2
pdf(file.path(output.dir, "heatmap_lncRNA.pdf"),width=5.5,height=5.5)
print(ht2)
dev.off()

##
library(MySeuratWrappers)  
seurat.obj <- readRDS("~/anno/harmony_anno_20_0.8.rds")
features <- c("PSMB1","ZFAS1","SNHG32","SNHG5","SNHG6","SNHG8","SNHG29","NEAT1","MALAT1","GAS5",
              "SNHG7","OTUD6B-AS1","PHF1","NORAD","RGS5","FTX","SNHG9","SNHG16","CYTOR",
              "EPB41L4A-AS1","LINC00963","PCAT19","MEG3","LINC02802","MIR99AHG",
              "MAGI2-AS3","MIR22HG","KCNQ1OT1","CARMN","LINC01615",
              "CRNDE","SERTAD4-AS1","MIR100HG","DANCR","CD27-AS1","PCED1B-AS1",
              "LINC01871","PRKCQ-AS1","LINC00623","LINC-PINT","LINC00924",
              "LINC01781","LINC00926","ILF3-DT","FAM215B","LINC02362","CHL1-AS2",
              "FOXD3-AS1","NSMCE1-DT") 
mycolor <- c("#00CD9B","#C3A5F7FF","#098B8B","#F8766D","#FFB8A4FF","#D8C6F3FF", 
             '#B53E2B',"#ACD45E",'#8C549C',"#9ECABE","#2F528F","#E3AD68")

vlo <- VlnPlot(seurat.obj, features = features,  
               stacked=T,
               pt.size=0,  
               cols = mycolor,             
               direction = "horizontal",      
               x.lab = '', y.lab = '')+ 
  theme(axis.text.x = element_blank(),
        axis.line = element_line(size = 0.1),
        axis.ticks.y = element_line(size = 0.1),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size =14,colour = "black"),
        axis.text.y = element_text(size=10,colour = "black"))  

pdf(file.path(output.dir, "vlnplot_lncRNA.pdf"),width=10,height=4)
print(vlo)
dev.off()

##
features <-c("MALAT1","SNHG5","PCAT19","PCED1B-AS1")

library(RColorBrewer)
p0 <- FeaturePlot(seurat.obj,features = "PCED1B-AS1",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 
p1 <- FeaturePlot(seurat.obj,features = "MALAT1",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 
p2 <- FeaturePlot(seurat.obj,features = "SNHG5",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 
p3 <- FeaturePlot(seurat.obj,features = "PCAT19",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 

ggsave(file.path(output.dir, "umap_PCED1B-AS1.pdf"),
       plot = p0,
       height = 2.5,
       width = 3)
ggsave(file.path(output.dir, "umap_MALAT1.pdf"),
       plot = p1,
       height = 2.5,
       width = 3)
ggsave(file.path(output.dir, "umap_SNHG5.pdf"),
       plot = p2,
       height = 2.5,
       width = 3)
ggsave(file.path(output.dir, "umap_PCAT19.pdf"),
       plot = p3,
       height = 2.5,
       width = 3)

##
features <-c("ACKR1","PLVAP","IL7R","TRAC")

library(RColorBrewer)
p0 <- FeaturePlot(seurat.obj,features = "ACKR1",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 
p1 <- FeaturePlot(seurat.obj,features = "PLVAP",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 
p2 <- FeaturePlot(seurat.obj,features = "IL7R",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 
p3 <- FeaturePlot(seurat.obj,features = "TRAC",reduction = "tsne",pt.size = 0.05)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        plot.title = element_text(hjust = 0.5,color="black",size=10,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) 

ggsave(file.path(output.dir, "umap_ACKR1.pdf"),
       plot = p0,
       height = 2.5,
       width = 3)
ggsave(file.path(output.dir, "umap_PLVAP.pdf"),
       plot = p1,
       height = 2.5,
       width = 3)
ggsave(file.path(output.dir, "umap_IL7R.pdf"),
       plot = p2,
       height = 2.5,
       width = 3)
ggsave(file.path(output.dir, "umap_TRAC.pdf"),
       plot = p3,
       height = 2.5,
       width = 3)

##
seurat.obj <- readRDS("./anno/harmony_anno_20_0.8.rds")
features <- read.csv("./results-V23-Figure/Fig1/lncRNA/4-final_lncRNA.csv")
features <- features$lncRNA

counts <- as.data.frame(seurat.obj[["RNA"]]@counts)

#
output <- vector("double", ncol(counts))   
counts$gene <- rownames(counts) 

for (i in 1:11756) {                       
  output[[i]] <-{
    a <- counts[,c(i,11757)] 
    a[a==0] <- NA            
    b <- na.omit(a)
    
    f <- t(b)      
    f <- length(f)
    f <- f/2
    
    c <- subset(b,subset = rownames(b) %in% features) 
    c <- t(c)   
    z = length(c)
    z = z/2
    
    w <- c(colnames(b),z,f)
    w <- as.data.frame(w)}
}

output <- as.data.frame(output)
output <- output[-c(2),]
rownames(output) <- c("ID","lncRNA","mRNA") 

##
output <- t(output) 
class(output)

output <- as.data.frame(output)
output$group <- "lncRNA" 
output$lncRNA <- as.numeric(output$lncRNA)  

p1 <- ggplot(output,aes(x=group, y=lncRNA)) + 
  geom_jitter(aes(color=group),width=0.08,size=0.1,color="grey")+     
  geom_boxplot(size = 0.4, 
               width = 0.3, 
               color = '#f8a6ac' ,  
               alpha = 0 ) +
  labs(x='',y='Number of lncRNAs \n(detection per cell)',title='')+mytheme  
p1
pdf(file.path(output.dir, "boxplot_lncRNA.pdf"),width=1.2,height=3)
print(p1)
dev.off()

##
output$group <- "mRNA" 
output$mRNA <- as.numeric(output$mRNA)  

p2 <- ggplot(output,aes(x=group, y=mRNA)) + 
  geom_jitter(aes(color=group),width=0.08,size=0.1,color="grey")+     
  geom_boxplot(size = 0.4, 
               width = 0.3, 
               color =  "#0072B2",  
               alpha = 0 )+ 
  labs(x='',y='Number of mRNAs \n(detection per cell)',title='') + mytheme 
p2
pdf(file.path(output.dir, "boxplot_mRNA.pdf"),width=1.4,height=3)
print(p2)
dev.off()

##
library(Seurat)
library(tidyverse)
library(ggplot2)
output.dir <- paste0("~/results-V23-Figure/Fig1") 
mytheme <- theme(
  axis.ticks=element_line(color="black",size=0.1,lineend = 10), 
  axis.title=element_text(color="black",size=10), 
  axis.text.y=element_text(size=8, color = "black"),  
  axis.text.x=element_text(size=10, color = "black"), 
  axis.line = element_line(linetype=1,color="black",size = 0.05),
  panel.grid = element_blank(), 
  legend.position = c('none'),  
  panel.background=element_rect(fill="white",colour="black",size=0.25)) 


##
output <- vector("double", ncol(counts))  
counts$gene <- rownames(counts) 

for (i in 1:11756) {                       
  output[[i]] <-{
    a <- counts[,c(i,11757)] 
    a[a==0] <- NA            
    b <- na.omit(a) 
    
    c1 <- t(b)
    d1 <- length(c1)
    d1 <- d1/2 
    
    c3 <- subset(b, subset = b[,1] > 3)  
    d3 <- t(c3)      
    e3 <- length(d3)
    e3 <- e3/2 
    f <- e3/d1
    
    c4 <- subset(b,subset = rownames(b) %in% features) 
    c5 <- subset(c4, subset = c4[,1] > 3)
    d5 <- t(c5)   
    e5 <- length(d5)
    e5 <- e5/2 
    # g <- e5/e4 ####counts＞3的lncRNA的比例
    g <- e5/d1
    
    w <- c(colnames(b),f,g)
    w <- as.data.frame(w)}
}

output <- as.data.frame(output)
output <- output[-c(2),] 
rownames(output) <- c("ID","mRNA","lncRNA")

#  
output <- t(output) 
class(output) 
output <- as.data.frame(output) 

output$lncRNA <- as.numeric(output$lncRNA)  
output$mRNA <- as.numeric(output$mRNA) 

#
output1 <- as.data.frame(output[,2])
output1$gene <- "mRNA"
colnames(output1) <- c("num","mRNA")
#
output2 <- as.data.frame(output[,3])
output2$gene <- "lncRNA"
colnames(output2) <- c("num","lncRNA")

##
colnames(output1) <- c("num","group")
colnames(output2) <- c("num","group")
A <- rbind(output1,output2)

#
library(RColorBrewer)
library(ggpubr)
A$group <- factor(A$group, levels=c("mRNA","lncRNA"))
p1 <- ggplot(A, aes(x = factor(group, levels = c("mRNA","lncRNA")), y = num, color = group)) + 
  geom_jitter(aes(color=group),width=0.08,size=0.1,color="grey")+  
  geom_boxplot(size = 0.4, 
               width = 0.3, 
               alpha = 0 )+
  stat_compare_means(label = "p.signif",comparisons = list(c("mRNA", "lncRNA")),method="wilcox.test")+#添加检验
  scale_color_manual(values=c( "#0072B2",'#f8a6ac'))+  
  ylim(0, 0.3)+
  labs(x='',y='Frequency of detection \n(gene >3 counts)',title='')+
  mytheme+
  NoLegend()  
p1
pdf(file.path(output.dir, "boxplot_mRNA-lncRNA-V2.pdf"),width=2,height=3)
print(p1)
dev.off()


##
p1 <- DimPlot(seurat_obj,
              reduction = "tsne",
              label = F,
              pt.size = 0.5,
              group.by = "cell_type",
              cols = color_cluster)+  
  ggtitle("mRNA")+
  theme(plot.title = element_text(hjust=0.5,size = 10), 
        legend.title = element_text(hjust=0.5),
        legend.text = element_text(size = 10), 
        axis.title = element_text(hjust=0.5,size = 10),
        axis.line = element_line(colour = "black",size = 0.01),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        element_line(colour = "black",size = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_color_manual(values=color_cluster)+
  guides(colour = guide_legend(override.aes = list(size=3), 
                               nrow = 12,
                               byrow = T,
                               reverse = F))  

pdf(file.path(output.dir,"tsne_anno_mRNA.pdf"),width=4.2,height=3)
print(p1)
dev.off()

##
gene_mata_mito <- read_csv("./results-V23-Figure/Fig1/lncRNA/real_exist_lncRNA.csv")
DefaultAssay(seurat_obj) <- "RNA" 
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
scale.genes <-  rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = scale.genes)
seurat_obj <- RunPCA(object = seurat_obj, pc.genes = gene_mata_mito)

ElbowPlot(seurat_obj, ndims = 50)
seurat_obj <- FindNeighbors(seurat_obj, 
                            reduction = "harmony", 
                            dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, 
                           resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, 
                      reduction = "harmony", 
                      dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, 
                      #reduction = "harmony", 
                      dims = 1:20)
library(ggsci)
set.seed(711)
p2 <- DimPlot(seurat_obj,
              reduction = "tsne",
              label = F,
              pt.size = 0.5,
              group.by = "cell_type",
              cols = color_cluster) + 
  ggtitle("lncRNA")+
  theme(plot.title = element_text(hjust=0.5,size = 10), 
        legend.title = element_text(hjust=0.5),
        legend.text = element_text(size = 10), 
        axis.title = element_text(hjust=0.5,size = 10),
        axis.line = element_line(colour = "black",size = 0.01),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),#加边框
        element_line(colour = "black",size = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_color_manual(values=color_cluster)+
  guides(colour = guide_legend(override.aes = list(size=3), 
                               nrow = 12,
                               byrow = T,
                               reverse = F))  

pdf(file.path(output.dir,"tsne_anno_lncRNA.pdf"),width=4.2,height=3)
print(p2)
dev.off()
