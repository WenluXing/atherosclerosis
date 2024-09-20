##
library(Seurat)
library(tidyverse)
library(ggplot2)

##
seurat.obj <- readRDS("./anno/harmony_anno_20_0.8.rds")
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
features <- read.csv("./results-V23-Figure/Fig1/lncRNA/4-final_lncRNA.csv")
features <- features$lncRNA

counts <- as.data.frame(seurat.obj[["RNA"]]@counts)

#for
output <- vector("double", ncol(counts))   
counts$gene <- rownames(counts)  

for (i in 1:11756) {                       
  output[[i]] <-{
    a <- counts[,c(i,11757)] #
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

#
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
pdf(file.path(output.dir, "boxplot_lncRNA.pdf"),width=1,height=2.2)
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
pdf(file.path(output.dir, "boxplot_mRNA.pdf"),width=1.15,height=2.2)
print(p2)
dev.off()

##
#for循环 
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
    g <- e5/d1
    
    w <- c(colnames(b),f,g)
    w <- as.data.frame(w)}
}

output <- as.data.frame(output)
output <- output[-c(2),] 
rownames(output) <- c("ID","mRNA","lncRNA") 

output <- t(output) 
class(output) 
output <- as.data.frame(output) 

output$lncRNA <- as.numeric(output$lncRNA)  
output$mRNA <- as.numeric(output$mRNA)  


output1 <- as.data.frame(output[,2])
output1$gene <- "mRNA"
colnames(output1) <- c("num","mRNA")

output2 <- as.data.frame(output[,3])
output2$gene <- "lncRNA"
colnames(output2) <- c("num","lncRNA")

#
colnames(output1) <- c("num","group")
colnames(output2) <- c("num","group")
A <- rbind(output1,output2)

library(RColorBrewer)
library(ggpubr)
A$group <- factor(A$group, levels=c("mRNA","lncRNA"))
p1 <- ggplot(A, aes(x = factor(group, levels = c("mRNA","lncRNA")), y = num, color = group)) + 
  geom_jitter(aes(color=group),width=0.08,size=0.1,color="grey")+  
  geom_boxplot(size = 0.4, 
               width = 0.3, 
               alpha = 0 )+
  stat_compare_means(label = "p.signif",comparisons = list(c("mRNA", "lncRNA")),method="wilcox.test")+
  scale_color_manual(values=c( "#0072B2",'#f8a6ac'))+ 
  ylim(0, 0.3)+ 
  labs(x='',y='Frequency of detection \n(gene >3 counts)',title='')+
  mytheme+
  NoLegend()  
p1
pdf(file.path(output.dir, "boxplot_mRNA-lncRNA.pdf"),width=1.2,height=2.4)
print(p1)
dev.off()


color_cluster <- c("#E64B35FF","#4DBBD5FF", "#00A087FF", 
                   "#3C5488FF", "#F39B7FFF", "#8491B4FF",
                   "#91D1C2FF", "#DC0000FF", "#7E6148FF", 
                   "#B09C85FF","#A2005699","#9467BDFF","#42B54099")

##
p1 <- DimPlot(seurat.obj,
              reduction = "tsne",
              label = F,
              pt.size = 0.5,
              group.by = "cell_type",
              cols = color_cluster)+  
  ggtitle("mRNA")+
  theme(plot.title = element_text(hjust=0.5,size = 8), 
        legend.title = element_text(hjust=0.5),
        legend.text = element_text(size = 8), 
        legend.position='right',
        axis.title = element_text(hjust=0.5,size = 8),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.25, linetype="solid"),
        element_line(colour = "black",size = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_color_manual(values=color_cluster)+
  guides(colour = guide_legend(override.aes = list(size=1), 
                               nrow = 12,
                               byrow = T,
                               reverse = F))  
p1
pdf(file.path(output.dir,"tsne_anno_mRNA.pdf"),width=2.8,height=2)
print(p1)
dev.off()

##
gene_mata_mito <- read_csv("./results-V23-Figure/Fig1/lncRNA/real_exist_lncRNA.csv")
DefaultAssay(seurat.obj) <- "RNA" 
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
scale.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = scale.genes)
seurat.obj <- RunPCA(object = seurat.obj, pc.genes = gene_mata_mito)

ElbowPlot(seurat.obj, ndims = 50)
seurat.obj <- FindNeighbors(seurat.obj,
                            reduction = "harmony",
                            dims = 1:20)
seurat.obj <- FindClusters(seurat.obj,
                           resolution = 0.8)
seurat.obj <- RunUMAP(seurat.obj,
                      reduction = "harmony",
                      dims = 1:20)
seurat.obj <- RunTSNE(seurat.obj,
                      dims = 1:20)
set.seed(711)
p2 <- DimPlot(seurat.obj,
              reduction = "tsne",
              label = F,
              pt.size = 0.5,
              group.by = "cell_type",
              cols = color_cluster) + 
  ggtitle("lncRNA")+
  theme(plot.title = element_text(hjust=0.5,size = 8), 
        legend.title = element_text(hjust=0.5),
        legend.text = element_text(size = 8), 
        legend.position='right',
        axis.title = element_text(hjust=0.5,size = 8),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.25, linetype="solid"),
        element_line(colour = "black",size = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_color_manual(values=color_cluster)+
  guides(colour = guide_legend(override.aes = list(size=1), 
                               nrow = 12,
                               byrow = T,
                               reverse = F))  
p2
pdf(file.path(output.dir,"tsne_anno_lncRNA.pdf"),width=2.8,height=2)
print(p2)
dev.off()


##
fTau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    print("No data for this gene avalable.")
  }
  return(res)
}

seurat.obj <- readRDS("./anno/harmony_anno_20_0.8.rds")
all_expr <- AverageExpression(seurat.obj,
                              assays = "RNA",
                              slot = "counts",
                              group.by = "cell_type")[["RNA"]]
FPKM <- as.data.frame(all_expr)

output <- c(rep(0,nrow(FPKM)))  
for (i in 1:nrow(FPKM)) {
  output[i] <- fTau(FPKM[i,])
}

output <- as.data.frame(output)
gene <- as.data.frame(rownames(FPKM))
Tau <- cbind(gene, output)
colnames(Tau) <- c("gene","Tau")

# 对Tau-gene分类
lncRNAs <- read.csv("./results-V23-Figure/Fig1/lncRNA/real_exist_lncRNA.csv")
Tau$group ="mRNA"

#提取
Tau1 <- Tau %>% filter(gene %in% (lncRNAs$lncRNA))
Tau1$group <- "lncRNA"
Tau2 <- rbind(Tau,Tau1)
table(Tau2$group)

p1 <- ggplot(Tau2, aes(x=factor(group, levels = c("mRNA","lncRNA")), y = Tau, fill = group)) +
  geom_violin(size=0.3) +
  geom_boxplot(width = 0.15, size=0.3, fill = "white", color = "black", alpha = 0) +
  scale_fill_manual(values = c('#f8a6ac','#95c2e7'))+
  labs(x = "", y = "Tau")+
  stat_compare_means(comparisons = list(c("mRNA","lncRNA")),
                     method = "wilcox.test", size=4,
                     aes(label = paste0("p : ", after_stat(p.signif)),),
                     label.y = 1) + mytheme +NoLegend()
p1
pdf(file.path(output.dir,"Vlnplot_Tau-han-count.pdf"),width=1.8,height=2.2)
print(p1)
dev.off()

##
library(MySeuratWrappers)  
features <- c("ZFAS1","SNHG5","SNHG6","SNHG8","SNHG29","NEAT1","MALAT1","GAS5",
              "SNHG7","OTUD6B-AS1","NORAD","FTX","SNHG9","SNHG16","CYTOR",
              "EPB41L4A-AS1",
              
              "LINC00963","PCAT19","MEG3","LINC02802","MIR99AHG",
              "MAGI2-AS3","MIR22HG","KCNQ1OT1","CARMN","LINC01615",
              "CRNDE","SERTAD4-AS1","MIR100HG","DANCR","CD27-AS1","PCED1B-AS1",
              "LINC01871","PRKCQ-AS1","LINC00623","LINC-PINT","LINC00924",
              "LINC01781","LINC00926","ILF3-DT","FAM215B","LINC02362","CHL1-AS2",
              "FOXD3-AS1","NSMCE1-DT") 
mycolor <- c("#E64B35FF","#4DBBD5FF","#00A087FF", 
             "#3C5488FF","#F39B7FFF","#8491B4FF",
             "#91D1C2FF","#DC0000FF","#7E6148FF", 
             "#B09C85FF","#A2005699","#9467BDFF","#42B54099")

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
        axis.title = element_text(size =8,colour = "black"),
        axis.text.y = element_text(size=8,colour = "black"))  
vlo
pdf(file.path(output.dir, "vlnplot_lncRNA.pdf"),width=8,height=3.1)
print(vlo)
dev.off()

features <-c("MALAT1","SNHG5","PCAT19","ACKR1","PCED1B-AS1","IL7R")
library(RColorBrewer)
p1 <- FeaturePlot(seurat.obj,features = "MALAT1",reduction = "tsne",pt.size = 0.03)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #去掉背景线
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,color="black",size=8,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) #加边框 #改变标题位置和字体大小
p2 <- FeaturePlot(seurat.obj,features = "SNHG5",reduction = "tsne",pt.size = 0.03)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #去掉背景线
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,color="black",size=8,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) #加边框 #改变标题位置和字体大小
p3 <- FeaturePlot(seurat.obj,features = "PCAT19",reduction = "tsne",pt.size = 0.03)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #去掉背景线
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,color="black",size=8,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) #加边框 #改变标题位置和字体大小

p4 <- FeaturePlot(seurat.obj,features = "ACKR1",reduction = "tsne",pt.size = 0.03)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #去掉背景线
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,color="black",size=8,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) #加边框 #改变标题位置和字体大小
p5 <- FeaturePlot(seurat.obj,features = "PCED1B-AS1",reduction = "tsne",pt.size = 0.03)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #去掉背景线
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,color="black",size=8,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) #加边框 #改变标题位置和字体大小
p6 <- FeaturePlot(seurat.obj,features = "IL7R",reduction = "tsne",pt.size = 0.03)+
  scale_color_gradient(low = "grey", high = "firebrick3")+
  scale_x_continuous("")+scale_y_continuous("")+
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #去掉背景线
        axis.ticks = element_blank(),
        axis.text = element_blank(), 
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,color="black",size=8,face = 'italic'),
        panel.border = element_rect(fill=NA,
                                    color="black", 
                                    size=0.2, 
                                    linetype="solid")) #加边框 #改变标题位置和字体大小

library(patchwork)
c <- p1+p2+p3+p4+p5+p6+plot_layout(ncol = 6)
ggsave(file.path(output.dir, "umap_lncRNA-mRNA.pdf"),
       plot = c,
       height = 2.2,
       width = 7.2)


##
source("~/share/20220802_scrna-m6A/custom_function.R")
source("~/share/20220802_scrna-m6A/custom_plot_function.R")

pathways <- c(
  "GOBP_CELL_CYCLE",
  "GOBP_REGULATION_OF_BMP_SIGNALING_PATHWAY",
  "GOBP_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN",
  "GOBP_VASCULOGENESIS",
  "GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_SIGNALING_PATHWAY",
  "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION",
  "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION",
  "GOBP_MONOCYTE_CHEMOTACTIC_PROTEIN_1_PRODUCTION",  
  "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
  "GOBP_CELL_MATRIX_ADHESION",
  "GOBP_LIPID_CATABOLIC_PROCESS",
  "GOBP_ENDOTHELIN_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
  "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
  "GOBP_INTERLEUKIN_2_PRODUCTION", 
  "GOBP_INTERLEUKIN_12_PRODUCTION", 
  "GOBP_INTERLEUKIN_23_PRODUCTION",
  "GOBP_INTERLEUKIN_6_MEDIATED_SIGNALING_PATHWAY", 
  "GOBP_RESPONSE_TO_INTERFERON_GAMMA"
)

End_GO <- read.csv("./Endothelial/k_20/correlation/GO:BP.cor.csv")
End_GO$cell_type <- "Endothelial"
End_selected_genes <- c("MALAT1","GAS5","SNHG29") 
df_End <- End_GO[End_GO$feature_y %in% pathways,]
df_End <- df_End[df_End$feature_x %in% End_selected_genes,]
SMC_GO <- read.csv("./SMC/k_20/correlation/GO:BP.cor.csv")
SMC_GO$cell_type <- "SMC"
SMC_selected_genes <- c("MALAT1","MIR100HG","NEAT1","SNHG8","SNHG9","SNHG29") 
df_SMC <- SMC_GO[SMC_GO$feature_y %in% pathways,]
df_SMC <- df_SMC[df_SMC$feature_x %in% SMC_selected_genes,]
F_GO <- read.csv("./Fibroblast/k_20/correlation/GO:BP.cor.csv")
F_GO$cell_type <- "Fibroblast"
F_selected_genes <- c("MALAT1","NEAT1","SNHG5","SNHG9","SNHG29","SERTAD4-AS1") #
df_F <- F_GO[F_GO$feature_y %in% pathways,]
df_F <- df_F[df_F$feature_x %in% F_selected_genes,]
Mac_selected_genes <- c("DANCR","GAS5","SNHG8","FTX","NEAT1","SENCR","MEXIS","LEXIS",
                        "MALAT1","SNHG9","OTUD6B-AS1") 
Mac_GO <- read.csv("./Macrophage/k_20/correlation/GO:BP.cor.csv")
Mac_GO$cell_type <- "Macrophage"
df_Mac <- Mac_GO[Mac_GO$feature_y %in% pathways,]
df_Mac <- df_Mac[df_Mac$feature_x %in% Mac_selected_genes,]

b <- rbind(df_End, df_SMC, df_F, df_Mac)
b <- b %>%
  filter(p_value < 0.05,
         feature_y %in% pathways)

b$feature_y <- b$feature_y |>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
b$feature_y <- factor(b$feature_y, 
                      levels=pathways|>
                        str_replace_all("GOBP_", "") |>
                        str_replace_all("_", " ") |>
                        str_to_sentence())
type1 <- c("Endothelial","SMC","Fibroblast","Macrophage")
b$cell_type <- factor(b$cell_type, levels=type1)

#
p <- b  %>% 
  ggplot(aes(x = feature_x, y = feature_y))+
  geom_point(aes(size = -log10(p_value), color = estimate), shape = 15)+
  facet_grid(. ~ cell_type,
             space = "free",
             scales = "free") +
  scale_color_gradient2(breaks=c(-0.5,0,0.5),
                        low = "navy", mid = "white", high = "firebrick3")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="black"),
        legend.position = "top",
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.spacing.x = unit(0, "pt"),
        #panel.grid = element_line(size = 0.2, color = "lightgrey"),
        axis.text.x = element_text(size = 9, color = "black",face = "italic",angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title = element_blank(),
        plot.title=element_text(size = 9, color = "black"))
p
ggsave(file.path(output.dir, "dotplot_lncRNA_way_cor.pdf"),
       plot = p,
       height = 5.5,
       width = 8.5)

##
A <- read_csv("number_T.csv")
B <- read_csv("number_NK.csv")
C <- read_csv("number_Pla.csv")
D <- read_csv("number_B.csv")
E <- read_csv("number_Mac.csv")

H <- read_csv("number_Fib.csv")
I <- read_csv("number_Fim.csv")
J <- read_csv("number_Per.csv")
K <- read_csv("number_Neu.csv")

L <- read_csv("number_End.csv")
M <- read_csv("number_SMC.csv")
N <- read_csv("number_Mast.csv")

output <- rbind( A, B, C, D, E, H, I, J, K, L, M, N)
colnames(output)[colnames(output) == "gene"] <- "mRNA"
output$cell_type <- factor(output$cell_type,
                           levels = c("Fibroblast","Fibromyocyte","SMC",
                                      "Endothelial","Macrophage","Plasma cell","Pericyte",
                                      "NK cell","Neuron","Mast cell","T cell","B cell")) 

#
p1 <- ggplot(output,aes(x=cell_type, y=mRNA)) + 
  geom_jitter(aes(color=group1),width=0.08,size=0.1,color="grey")+  
  geom_boxplot(size = 0.4, 
               width = 0.4, 
               color = '#95c2e7',
               alpha = 0)+
  ylim(0, 4000)+ 
  labs(x='',y='number of mRNAs',title='')+ 
  theme_classic() + mytheme + 
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1)) 

pdf(file.path(output.dir,"FigS1B_boxplot_celltype_mRNA.pdf"),width=2.8,height=2.5)
print(p1)
dev.off()

#
p2 <- ggplot(output,aes(x=cell_type, y=lncRNA)) + 
  geom_jitter(aes(color=group2),width=0.08,size=0.1,color="grey")+   
  geom_boxplot(size = 0.4, 
               width = 0.4, 
               color = '#f8a6ac', 
               alpha = 0)+
  ylim(0, 100)+ 
  labs(x='',y='number of lncRNAs',title='')+
  theme_classic() + mytheme + 
  theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1)) 

pdf(file.path(output.dir,"boxplot_celltype_lncRNA.pdf"),width=2.8,height=2.5)
print(p2)
dev.off()

library(circlize)
colors <- colorRampPalette(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50) ##蓝到红
values <- seq(-1, 1, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)

##
Bcell <- readRDS("./B cell/k_20/data/seurat_obj.rds")
Bcell1 <- hdWGCNA::GetMetacellObject(Bcell)
End <- readRDS("./Endothelial/k_20/data/seurat_obj.rds")
End1 <- hdWGCNA::GetMetacellObject(End)
Fib <- readRDS("./Fibroblast/k_20/data/seurat_obj.rds")
Fib1 <- hdWGCNA::GetMetacellObject(Fib)
Fibm <- readRDS("./Fibromyocyte/k_20/data/seurat_obj.rds")
Fibm1 <- hdWGCNA::GetMetacellObject(Fibm)
Mac <- readRDS("./Macrophage/k_20/data/seurat_obj.rds")
Mac1 <- hdWGCNA::GetMetacellObject(Mac)
NK <- readRDS("./NK cell/k_20/data/seurat_obj.rds")
NK1 <- hdWGCNA::GetMetacellObject(NK)
Pericyte <- readRDS("./Pericyte/k_20/data/seurat_obj.rds")
Pericyte1 <- hdWGCNA::GetMetacellObject(Pericyte)
Plasma <- readRDS("./Plasma cell/k_20/data/seurat_obj.rds")
Plasma1 <- hdWGCNA::GetMetacellObject(Plasma)
SMC <- readRDS("./SMC/k_20/data/seurat_obj.rds")
SMC1 <- hdWGCNA::GetMetacellObject(SMC)
Tcell <- readRDS("./T cell/k_20/data/seurat_obj.rds")
Tcell1 <- hdWGCNA::GetMetacellObject(Tcell)

Neu <- readRDS("./Neuron/k_20/data/seurat_obj.rds")
Neu1 <- hdWGCNA::GetMetacellObject(Neu)

sce <- merge(Bcell1,c(End1,Fib1,Fibm1,Mac1,NK1,Pericyte1,Plasma1,SMC1,Tcell1,Neu1))
sce_matrix <- AverageExpression(sce,group.by = "donor",assays = "RNA",slot = "data")
sce_matrix <- sce_matrix$RNA

library(clusterProfiler)
CHRONIC_INFLAMMATORY_RESPONSE <- read.gmt("./GOBP_CHRONIC_INFLAMMATORY_RESPONSE.v2023.1.Hs.gmt")

sce_matrix_cor <- rcorr(t(sce_matrix[rownames(sce_matrix) %in% c(features,
                                                                 CHRONIC_INFLAMMATORY_RESPONSE$gene),]))
sce_matrix_cor_r <- sce_matrix_cor$r
sce_matrix_cor_p <- sce_matrix_cor$P

lncRNA_cor <- sce_matrix_cor_r[rownames(sce_matrix_cor_r) %in% CHRONIC_INFLAMMATORY_RESPONSE$gene,features]
lncRNA_cor[is.na(lncRNA_cor)] <- 0 
lncRNA_p <- sce_matrix_cor_p[rownames(lncRNA_cor),colnames(lncRNA_cor)]

library(stringr)
library(pheatmap)

lncRNA_cor <- t(lncRNA_cor)
ht1 <- pheatmap(lncRNA_cor,
                cluster_cols = T,treeheight_col = 0,
                cluster_rows = T,treeheight_row = 15,
                show_rownames = T, 
                col = colors ,
                angle_col = 90, 
                main="CHRONIC_INFLAMMATORY_RESPONSE"|>
                  str_replace_all("_", " ") |>
                  str_to_sentence(),
                fontsize = 9) 
ht1
gene <- c("NEAT1","MALAT1","GAS5","LINC02802","MIR99AHG","PCED1B-AS1","MEG3","KCNQ1OT1","ZFAS1",
          "LINC01781","LINC-PINT","LINC02362","SNHG9","MEXIS","SCNER")
pdf(file.path(output.dir,"FigS1C_Chronic_lncRNA_cor.pdf"),height = 10, width = 5)
ht2 <- add.flag(ht1, kept.labels = gene, repel.degree = 0.2) #！！！
print(ht2)
dev.off() 

ht2$grobs[[4]]$gp=grid::gpar(fontface="italic")
ht2$grobs[[5]]$gp=grid::gpar(fontface="italic")
plot(ht2)

pdf(file.path(output.dir,"Chronic_lncRNA_cor.pdf"),height = 10, width = 5)
grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
grid.draw(ht2)
dev.off() 

##
lncRNA <- c("SNHG29","GAS5","SNHG8","NEAT1","MALAT1",
            "SNHG9","CYTOR","MEG3","KCNQ1OT1","CARMN","SERTAD4-AS1","MIR100HG",
            "DANCR","ILF3-DT","ZFAS1","SNHG5")

End_GO <- read.csv("./Endothelial/k_20/correlation/gene.cor.csv")
End_GO$cell_type <- "Endothelial"
End_selected_genes <- c("ACKR1","CLDN5","PECAM1","ACE","EDNRB","PDGFD","PDGFRA","CD36")
df_End <- End_GO[End_GO$feature_y %in% End_selected_genes,]
df_End <- df_End[df_End$feature_x %in% lncRNA,]

SMC_GO <- read.csv("./SMC/k_20/correlation/gene.cor.csv")
SMC_GO$cell_type <- "SMC"
SMC_selected_genes <- c("ACTA2","MYH11","ACTG2","ACTN1","BACH1","BACH2","TCF21","PRMT5","IRF8","NOTCH3","KLF4","MYL9","MGP","VIM") 
df_SMC <- SMC_GO[SMC_GO$feature_y %in% SMC_selected_genes,]
df_SMC <- df_SMC[df_SMC$feature_x %in% lncRNA,]

F_GO <- read.csv("./Fibroblast/k_20/correlation/gene.cor.csv")
F_GO$cell_type <- "Fibroblast"
F_selected_genes <- c("DCN","LUM","MEOX1") 
df_F <- F_GO[F_GO$feature_y %in% F_selected_genes,]
df_F <- df_F[df_F$feature_x %in% lncRNA,]

Mac_selected_genes <- c("PLTP","APOE", "CCL3","CXCL16","CD83","CD86","TNF","HLA-A","IFNGR1","SLC11A1","JUN","JUNB"
) 
Mac_GO <- read.csv("./Macrophage/k_20/correlation/gene.cor.csv")
Mac_GO$cell_type <- "Macrophage"
df_Mac <- Mac_GO[Mac_GO$feature_y %in% Mac_selected_genes,]
df_Mac <- df_Mac[df_Mac$feature_x %in% lncRNA,]

b <- rbind(df_End, df_SMC, df_F, df_Mac)
b <- b %>%
  filter(p_value < 0.05,
         feature_x %in% lncRNA)

b$feature_x <- factor(b$feature_x,
                      levels=lncRNA)
b$cell_type <- factor(b$cell_type,
                      levels=c("Endothelial","SMC","Fibroblast","Macrophage"))

b$feature_y <- factor(b$feature_y,
                      levels=c(End_selected_genes,SMC_selected_genes,
                               F_selected_genes,Mac_selected_genes))
#
p <- b  %>%
  ggplot(aes(x = feature_y, y = feature_x))+
  geom_point(aes(size = -log10(p_value), color = estimate), shape = 10)+
  facet_grid(. ~ cell_type,
             space = "free",
             scales = "free") +
  scale_color_gradient2(#name="cor",
    breaks=c(-0.5,0,0.5),
    low = "navy", mid = "white", high = "firebrick3")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="black"),
        legend.position = "right",
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.spacing.x = unit(0, "pt"),
        #panel.grid = element_line(size = 0.2, color = "lightgrey"),
        axis.text.x = element_text(size = 8, color = "black",face = "italic",angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 8, color = "black",face = "italic"),
        axis.title = element_blank(),
        plot.title=element_text(size = 8, color = "black"))
p
ggsave(file.path(output.dir, "dotplot_lncRNA-marker_cor.pdf"),
       plot = p,
       height = 3.5,
       width = 8.8)

##
lncRNA <- read_csv("./results-V23-Figure/Fig1/lncRNA/real_exist_lncRNA.csv")
lncRNA <- lncRNA$lncRNA
adata_1 <- read_csv("./End/tfs_targer.tsv")
adata_1 <- adata_1 %>% filter(target_gene %in% lncRNA) %>% mutate(cell_type="Endothelial")
adata_2 <- read_csv("./SMC/tfs_targer.tsv")
adata_2 <- adata_2 %>% filter(target_gene %in% lncRNA)%>% mutate(cell_type="SMC")
adata_3 <- read_csv("./Mac/tfs_targer.tsv")
adata_3 <- adata_3 %>% filter(target_gene %in% lncRNA)%>% mutate(cell_type="Macrophage")
adata_4 <- read_csv("./Fib/tfs_targer.tsv")
adata_4 <- adata_4 %>% filter(target_gene %in% lncRNA)%>% mutate(cell_type="Fibroblast")
adata <- rbind(adata_1,adata_2,adata_3,adata_4)

scenic <- c("KLF2","KLF4","KLF6","KLF14","EGR1", 
            "GATA2","TAF1","HIF1A","JUNB","OCT4","NFKB1","NFKB2","ETS2","NR2F1","IRF1")
lncRNA <- c("MALAT1","NEAT1","MIR100HG","PCAT19","SNHG8","SNHG9","FTX",
            "XIST","SNHG5","SNHG6","SNHG16","SNHG29","GAS5","MEG3","KCNQ1OT1",
            "H19","ZFAS1","NORAD","CYTOR","DANCR","SENCR","CARMN","LUCAT1",
            "CDKN2B-AS1","EPB41L4A-AS1","SERTAD4-AS1","NEXN-AS1","OIP5-AS1","CEBPB-AS1",
            "PVT1","LINC-PINT","LINC00623","LINC00863","ADAMTS9-AS2","C2orf27A")
data <- adata %>%
  filter(target_gene %in% lncRNA) %>%
  filter(tf %in% scenic) %>%
  mutate(num = 10, a = 1) 

data$cell_type <- factor(data$cell_type, levels = c("Endothelial","SMC",'Fibroblast',"Macrophage"))
data$tf <- factor(data$tf, levels = scenic)
# 
p <- ggplot(data, aes(target_gene, tf))+
  geom_point(aes(fill=cell_type),color="black",size=3,shape=21,alpha=0.9,stroke=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, face = "italic", angle = 90, vjust=0.5, hjust=1, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.line = element_line(size = 0.1, color = "black"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.title = element_blank(),
        plot.title=element_text(size = 8, color = "black"),
        panel.border = element_rect(color = "black", size = 0.3, fill = NA),
        panel.grid =element_blank(),
        panel.spacing.x = unit(0, "pt"))+  
  facet_grid(. ~cell_type,
             space = "free",
             scales = "free") +
  NoLegend()
p
ggsave(file.path(output.dir,"tfs.pdf"),
       plot = p,
       height = 2.8,
       width = 8)
