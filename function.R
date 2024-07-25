#### 加载包 ----
library(Hmisc)
library(Seurat)
source("~/share/custom_function.R")
source("~/share/custom_plot_function.R")
output.dir <- paste0("~/results-V23-Figure/Fig2") 
library(circlize)
colors <- colorRampPalette(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(50) ##蓝到红
values <- seq(-1, 1, length.out = 51)[-51]
col_fun = colorRamp2(values, colors)

##
Bcell <- readRDS("~/results-V3-metacell/B cell/k_20/data/seurat_obj.rds")
Bcell1 <- hdWGCNA::GetMetacellObject(Bcell)
End <- readRDS("~/results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")
End1 <- hdWGCNA::GetMetacellObject(End)
Fib <- readRDS("~/results-V3-metacell/Fibroblast/k_20/data/seurat_obj.rds")
Fib1 <- hdWGCNA::GetMetacellObject(Fib)
Fibm <- readRDS("~/results-V3-metacell/Fibromyocyte/k_20/data/seurat_obj.rds")
Fibm1 <- hdWGCNA::GetMetacellObject(Fibm)
Mac <- readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")
Mac1 <- hdWGCNA::GetMetacellObject(Mac)
NK <- readRDS("~/results-V3-metacell/NK cell/k_20/data/seurat_obj.rds")
NK1 <- hdWGCNA::GetMetacellObject(NK)
Pericyte <- readRDS("~/results-V3-metacell/Pericyte/k_20/data/seurat_obj.rds")
Pericyte1 <- hdWGCNA::GetMetacellObject(Pericyte)
Plasma <- readRDS("~/results-V3-metacell/Plasma cell/k_20/data/seurat_obj.rds")
Plasma1 <- hdWGCNA::GetMetacellObject(Plasma)
SMC <- readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
SMC1 <- hdWGCNA::GetMetacellObject(SMC)
Tcell <- readRDS("~/results-V3-metacell/T cell/k_20/data/seurat_obj.rds")
Tcell1 <- hdWGCNA::GetMetacellObject(Tcell)
Neu <- readRDS("~/results-V3-metacell/Neuron/k_20/data/seurat_obj.rds")
Neu1 <- hdWGCNA::GetMetacellObject(Neu)

sce <- merge(Bcell1,c(End1,Fib1,Fibm1,Mac1,NK1,Pericyte1,Plasma1,SMC1,Tcell1,Neu1))
sce_matrix <- AverageExpression(sce,group.by = "donor",assays = "RNA",slot = "data")
sce_matrix <- sce_matrix$RNA

#
features <- read.csv("./results-V23-Figure/Fig1/lncRNA/real_exist_lncRNA.csv")
features <- features$lncRNA

#
library(clusterProfiler)
CHRONIC_INFLAMMATORY_RESPONSE <- read.gmt("~/results-V6-cor-lncRNA和way-Fig2/define_pathway_self/GOBP_CHRONIC_INFLAMMATORY_RESPONSE.v2023.1.Hs.gmt")

##
sce_matrix_cor <- rcorr(t(sce_matrix[rownames(sce_matrix) %in% c(features,
                                                                 CHRONIC_INFLAMMATORY_RESPONSE$gene),]))
sce_matrix_cor_r <- sce_matrix_cor$r
sce_matrix_cor_p <- sce_matrix_cor$P

lncRNA_cor <- sce_matrix_cor_r[rownames(sce_matrix_cor_r) %in% CHRONIC_INFLAMMATORY_RESPONSE$gene,features]
lncRNA_cor[is.na(lncRNA_cor)] <- 0 
lncRNA_p <- sce_matrix_cor_p[rownames(lncRNA_cor),colnames(lncRNA_cor)]

# 
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
                # fontsize_row = 8,  
                # fontsize_col = 10, 
                fontsize = 9) 
ht1
gene <- c("NEAT1","MALAT1","GAS5","LINC02802","MIR99AHG","PCED1B-AS1","MEG3","KCNQ1OT1","ZFAS1",
          "LINC01781","LINC-PINT","LINC02362","SNHG9","MEXIS","SCNER")
pdf(file.path(output.dir,"Chronic_lncRNA_cor.pdf"),height = 10,width = 5)
ht2 <- add.flag(ht1, kept.labels = gene, repel.degree = 0.2) #！！！
print(ht2)
dev.off() 

ht2$grobs[[4]]$gp=grid::gpar(fontface="italic")
ht2$grobs[[5]]$gp=grid::gpar(fontface="italic")
plot(ht2)

pdf(file.path(output.dir,"Chronic_lncRNA_cor.pdf"),height = 10,width = 5)
grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
grid.draw(ht2)
dev.off() 


##
INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS <- read.gmt("~/results-V6-cor-lncRNA和way-Fig2/define_pathway_self/GOBP_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS.v2023.1.Hs.gmt")
CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE <- read.gmt("~/results-V6-cor-lncRNA和way-Fig2/define_pathway_self/GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2023.1.Hs.gmt")

##
library(stringr)
library(pheatmap)

# 
sce_matrix_cor <- rcorr(t(sce_matrix[rownames(sce_matrix) %in% c(features,
                                                                 INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS$gene),]))
sce_matrix_cor_r <- sce_matrix_cor$r
sce_matrix_cor_p <- sce_matrix_cor$P

lncRNA_cor <- sce_matrix_cor_r[rownames(sce_matrix_cor_r) %in% INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS$gene,features]
lncRNA_cor[is.na(lncRNA_cor)] <- 0 
lncRNA_p <- sce_matrix_cor_p[rownames(lncRNA_cor),colnames(lncRNA_cor)]

ht2 <- pheatmap(lncRNA_cor,
                cluster_cols = T,treeheight_col = 0,
                cluster_rows = T,treeheight_row = 15,
                show_rownames = T, 
                show_colnames = F,
                color =colorRampPalette(colors=c("#010076","white","#D21B23"))(100),
                angle_col = 90,
                main="INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS"|>
                  str_replace_all("_", " ") |>
                  str_to_sentence(),
                fontsize = 9) 
gene <- c("NOTCH2","IL10","CXCR2","HMGB2","TNF","MAPK14","NOTCH1","GATA3","CCR7","TREM2")
pdf(file.path(output.dir,"Antigenic_lncRNA_cor.pdf"),height = 4,width = 10)
ht2 <- add.flag(ht2, kept.labels = gene, repel.degree = 0.2) 
print(ht2)
dev.off() 

ht2$grobs[[4]]$gp=grid::gpar(fontface="italic")
plot(ht2)

pdf(file.path(output.dir,"Antigenic_lncRNA_cor.pdf"),height = 4,width = 10)
grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
grid.draw(ht2)
dev.off() 

# 
sce_matrix_cor <- rcorr(t(sce_matrix[rownames(sce_matrix) %in% c(features,
                                                                 CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE$gene),]))
sce_matrix_cor_r <- sce_matrix_cor$r
sce_matrix_cor_p <- sce_matrix_cor$P

lncRNA_cor <- sce_matrix_cor_r[rownames(sce_matrix_cor_r) %in% CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE$gene,features]
lncRNA_cor[is.na(lncRNA_cor)] <- 0 
lncRNA_p <- sce_matrix_cor_p[rownames(lncRNA_cor),colnames(lncRNA_cor)]

ht3 <- pheatmap(lncRNA_cor,
                cluster_cols = T,treeheight_col = 0,
                cluster_rows = T,treeheight_row = 15,
                show_rownames = T, 
                show_colnames = F,
                color =colorRampPalette(colors=c("#010076","white","#D21B23"))(100),
                angle_col = 90,
                main="CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE"|>
                  str_replace_all("_", " ") |>
                  str_to_sentence(),
                fontsize = 9) 
gene <- c("IL1R2","APOD","IL17B","TNF","IL6","HIF1A","STAT3","MAPK14","TREM2","TLR4")
pdf(file.path(output.dir,"Cytokine_lncRNA_cor.pdf"),height = 4,width = 10)
ht3 <- add.flag(ht3, kept.labels = gene, repel.degree = 0.2) 
print(ht3)
dev.off() 

ht3$grobs[[4]]$gp=grid::gpar(fontface="italic")
plot(ht3)

pdf(file.path(output.dir,"Cytokine_lncRNA_cor.pdf"),height = 4,width = 10)
grid.draw(rectGrob(gp=gpar(fill="white", lwd=0)))
grid.draw(ht3)
dev.off() 


##
pathways <- c(
  "GOBP_CELL_CYCLE",
  "GOBP_REGULATION_OF_BMP_SIGNALING_PATHWAY",
  "GOBP_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN",
  "GOBP_VASCULOGENESIS",
  "GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_SIGNALING_PATHWAY",
  "GOBP_BLOOD_VESSEL_REMODELING",
  "GOBP_MONOCYTE_CHEMOTACTIC_PROTEIN_1_PRODUCTION",  
  "GOBP_PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_BETA_SIGNALING_PATHWAY",
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

End_GO <- read.csv("~/results-V4-correlation/Endothelial/k_20/correlation/GO:BP.cor.csv")
End_GO$cell_type <- "Endothelial"
End_selected_genes <- c("SNHG29","GAS5","MALAT1","MAIT","ANRIL") 
df_End <- End_GO[End_GO$feature_y %in% pathways,]
df_End <- df_End[df_End$feature_x %in% End_selected_genes,]

SMC_GO <- read.csv("~/results-V4-correlation/SMC/k_20/correlation/GO:BP.cor.csv")
SMC_GO$cell_type <- "SMC"
SMC_selected_genes <- c("SNHG32","SNHG29","NEAT1","MIR100HG","SNHG9","TUG1","ANRIL","SENCR",
                        "SNHG8","LINC02362","CHL1-AS2") 
df_SMC <- SMC_GO[SMC_GO$feature_y %in% pathways,]
df_SMC <- df_SMC[df_SMC$feature_x %in% SMC_selected_genes,]

F_GO <- read.csv("~/results-V4-correlation/Fibroblast/k_20/correlation/GO:BP.cor.csv")
F_GO$cell_type <- "Fibroblast"
F_selected_genes <- c("SNHG29","SNHG5","SNHG9","NEAT1","SERTAD4-AS1") #
df_F <- F_GO[F_GO$feature_y %in% pathways,]
df_F <- df_F[df_F$feature_x %in% F_selected_genes,]

Mac_selected_genes <- c("DANCR","GAS5","SNHG8","FTX","NEAT1","SENCR","MEXIS","LEXIS",
                        "MALAT1","SNHG9","OTUD6B-AS1") 
Mac_GO <- read.csv("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.cor.csv")
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

##
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

        axis.text.x = element_text(size = 9, color = "black",face = "italic",angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title = element_blank(),
        plot.title=element_text(size = 9, color = "black"))
p
ggsave(file.path(output.dir, "otplot_lncRNA_way_cor.pdf"),
       plot = p,
       height = 5.5,
       width = 8.5)


##
lncRNA <- c("ZFAS1","SNHG32","SNHG5","GAS5","SNHG29","SNHG8","NEAT1","MALAT1",
            "SNHG9","CYTOR","MEG3","KCNQ1OT1","CARMN","SERTAD4-AS1","MIR100HG",
            "DANCR","ILF3-DT","FAM215B","FOXD3-AS1")

End_GO <- read.csv("~/results-V4-correlation/Endothelial/k_20/correlation/gene.cor.csv")
End_GO$cell_type <- "Endothelial"
End_selected_genes <- c("ACKR1","CLDN5","VWF","PECAM1","ACE","EDNRB","PDGFD","PDGFRA","CD36") 
df_End <- End_GO[End_GO$feature_y %in% End_selected_genes,]
df_End <- df_End[df_End$feature_x %in% lncRNA,]

SMC_GO <- read.csv("~/results-V4-correlation/SMC/k_20/correlation/gene.cor.csv")
SMC_GO$cell_type <- "SMC"
SMC_selected_genes <- c("ACTA2","MYH11","ACTG2","ACTN1","BACH1","BACH2","TCF21","PRMT5","IRF8","NOTCH3","KLF4","MYL9","MGP","VIM") 
df_SMC <- SMC_GO[SMC_GO$feature_y %in% SMC_selected_genes,]
df_SMC <- df_SMC[df_SMC$feature_x %in% lncRNA,]

F_GO <- read.csv("~/results-V4-correlation/Fibroblast/k_20/correlation/gene.cor.csv")
F_GO$cell_type <- "Fibroblast"
F_selected_genes <- c("DCN","LUM","MEOX1") 
df_F <- F_GO[F_GO$feature_y %in% F_selected_genes,]
df_F <- df_F[df_F$feature_x %in% lncRNA,]

Mac_selected_genes <- c("C1QA","C1QB","C1QC","MIF","CD63","LGMN","PLTP") 
Mac_GO <- read.csv("~/results-V4-correlation/Macrophage/k_20/correlation/gene.cor.csv")
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
        legend.position = "top",
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.spacing.x = unit(0, "pt"),
        axis.text.x = element_text(size = 9, color = "black",face = "italic",angle = 90,hjust = 1,vjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black",face = "italic"),
        axis.title = element_blank(),
        plot.title=element_text(size = 9, color = "black"))
p
ggsave(file.path(output.dir, "dotplot_lncRNA-marker_cor.pdf"),
       plot = p,
       height = 5,
       width = 8)


##
lncRNA <- read_csv("./results-V23-Figure/Fig1/lncRNA/real_exist_lncRNA.csv")
lncRNA <- lncRNA$lncRNA

adata_1 <- read_csv("./results-V7-scenic/End/tfs_targer.tsv")
adata_1 <- adata_1 %>% filter(target_gene %in% lncRNA) %>% mutate(cell_type="Endothelial")

adata_2 <- read_csv("./results-V7-scenic/SMC/tfs_targer.tsv")
adata_2 <- adata_2 %>% filter(target_gene %in% lncRNA)%>% mutate(cell_type="SMC")

adata_3 <- read_csv("./results-V7-scenic/Mac/tfs_targer.tsv")
adata_3 <- adata_3 %>% filter(target_gene %in% lncRNA)%>% mutate(cell_type="Macrophage")

adata_4 <- read_csv("./results-V7-scenic/Fib/tfs_targer.tsv")
adata_4 <- adata_4 %>% filter(target_gene %in% lncRNA)%>% mutate(cell_type="Fibroblast")

adata <- rbind(adata_1,adata_2,adata_3,adata_4)
          
scenic <- c("KLF2","KLF3","KLF4","KLF6", "KLF14","GATA2",
            "EGR1","TAF1","SPI1","JUNB","OCT4", "NFKB1","NFKB2","ETS2", 
            "NR2F1","IRF1","FOS","HIF1A" )
lncRNA <- c("MALAT1","NEAT1","MIR100HG","PCAT19","SNHG8","SNHG9","FTX",
            "XIST","SNHG5","SNHG6","SNHG16","SNHG32","SNHG29","GAS5","MEG3","KCNQ1OT1",
            "H19","ZFAS1","PHF1","NORAD","CYTOR","DANCR","SENCR","MEXIS","CARMN","LUCAT1",
            "CDKN2B-AS1","EPB41L4A-AS1","SERTAD4-AS1","NEXN-AS1","OIP5-AS1","CEBPB-AS1","NEXN-AS1",
            "PVT1","LINC-PINT","LINC00623","LINC00863","ADAMTS9-AS2","C2orf27A")

data <- adata %>%
  filter(target_gene %in% lncRNA) %>%
  filter(tf %in% scenic)
data$cell_type <- factor(data$cell_type, levels = c("Endothelial","SMC",'Fibroblast',"Macrophage"))
data$tf <- factor(data$tf, levels = scenic)

p <- ggplot(data, aes(target_gene, tf))+
  geom_point(aes(fill=cell_type),color="black",size=4,shape=21,alpha=0.9,stroke=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, face = "italic", angle = 90, vjust=0.5, hjust=1, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.line = element_line(size = 0.1, color = "black"),
        axis.ticks = element_line(size = 0.3, color = "black"),
        axis.title = element_blank(),
        plot.title=element_text(size = 8, color = "black"),
        panel.border = element_rect(color = "black", size = 0.3, fill = NA),
        panel.grid =element_blank(),
        panel.spacing.x = unit(0, "pt"))+  #修改分面面板间距
  facet_grid(. ~cell_type,
             space = "free",
             scales = "free") +
  NoLegend()
p
ggsave(file.path(output.dir,"tfs.pdf"),
  plot = p,
  height = 3.6,
  width = 9.5)
