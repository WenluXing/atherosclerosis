#### 加载包 ----
library(Seurat)
source("~/share/custom_function.R")
source("~/share/custom_plot_function.R")
output.dir <- paste0("~/results-V24-Figure/Fig3") 
selected_cell_type <- "Endothelial"

##
selected_genes <- c("MALAT1","PCAT19","CYTOR","SERTAD4-AS1","GAS5","SNHG29","ZFAS1","SNHG5")
adata <- read.csv("~/results-V4-correlation/Endothelial/k_20/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value < 0.05)

selected_features <- c(
  "GOBP_ENDOTHELIAL_CELL_MIGRATION",
  "GOBP_ENDOTHELIAL_CELL_MORPHOGENESIS",
  "GOBP_ENDOTHELIAL_CELL_MATRIX_ADHESION",
  "GOBP_ENDOTHELIAL_CELL_PROLIFERATION",
  "GOBP_ANGIOTENSIN_ACTIVATED_SIGNALING_PATHWAY",
  "GOBP_VASCULOGENESIS",
  "GOBP_LEUKOCYTE_MIGRATION",  
  "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION",  
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_POSITIVE_CHEMOTAXIS",
  "GOBP_MACROPHAGE_CHEMOTAXIS",
  "GOBP_POSITIVE_REGULATION_OF_MONOCYTE_CHEMOTAXIS",
  "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION")
a <- filtered_adata %>% 
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features)
a$feature_x <- factor(a$feature_x, levels=selected_genes)
a$feature_y <- factor(a$feature_y, levels=selected_features)
a$feature_y <- a$feature_y |> 
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

p1 <- a  %>%
  ggplot(aes(x = feature_x, y = feature_y))+
  geom_point(aes(size = -log10(p_value), color = estimate), shape = 19)+
  scale_color_gradient2(name="estimate",breaks=c(-0.4,0,0.5),low="navy",mid="white",high="firebrick3")+
  scale_y_discrete(position="right")+
  theme_bw() +
  coord_equal()+
  theme(legend.key.size = unit(6, "pt"),
        panel.grid.major = element_blank(),
        axis.line = element_line(color = "black",size=0.1), 
        legend.title = element_text(size = 7, color = "black"),
        legend.text = element_text(size = 7, color = "black"),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8, color = "black",face = "italic",angle = 45,vjust=1, hjust=1),
        axis.text.y = element_text(size = 8, color = "black"))
p1
ggsave(file.path(output.dir, "dotplot_pathway_lnc.pdf"),
       plot = p1,
       height = 3.5,
       width = 5.5)

##
selected_genes <- c("MALAT1","PCAT19","CYTOR","SERTAD4-AS1","GAS5","SNHG29","ZFAS1","SNHG5")
seurat_obj <- readRDS("~/results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")

table(seurat_obj$cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <- readRDS("~/results-V4-correlation/Endothelial/k_20/correlation/GO:BP.scores.rds")
scores[1:2, 1:2]

if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]

#
library(ggExtra)
P <- expr %>%
  ggplot(aes(x = MALAT1, y = GOBP_POSITIVE_CHEMOTAXIS)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MALAT1, y = GOBP_POSITIVE_CHEMOTAXIS),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_POSITIVE_CHEMOTAXIS" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "orange"),
                 yparams = list(fill = "#00AFBB"))
ggsave(file.path(output.dir, "san_MALAT1_chemotaxis.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = PCAT19, y = GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = PCAT19, y = GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_\nCELL_CELL_ADHESION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "orange"),
                 yparams = list(fill = "#00AFBB"))
ggsave(file.path(output.dir, "san_PCAT19_adhension.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = MALAT1, y = GOBP_VASCULOGENESIS)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MALAT1, y = GOBP_VASCULOGENESIS),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_VASCULOGENESIS" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "orange"),
                 yparams = list(fill = "#00AFBB"))
ggsave(file.path(output.dir, "san_MALAT1_vasculogenesis.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = PCAT19, y = GOBP_LEUKOCYTE_CELL_CELL_ADHESION)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = PCAT19, y = GOBP_LEUKOCYTE_CELL_CELL_ADHESION),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_LEUKOCYTE_CELL_CELL_ADHESION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "orange"),
                 yparams = list(fill = "#00AFBB"))
ggsave(file.path(output.dir, "san_PCAT19_adhension2.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = MALAT1, y = GOBP_POSITIVE_REGULATION_OF_MONOCYTE_CHEMOTAXIS)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MALAT1, y = GOBP_POSITIVE_REGULATION_OF_MONOCYTE_CHEMOTAXIS),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_POSITIVE_REGULATION_OF_\nMONOCYTE_CHEMOTAXIS" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "orange"),
                 yparams = list(fill = "#00AFBB"))
ggsave(file.path(output.dir, "san_MALAT1_chemotaxis2.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)

##
adata <-
  read_csv("~/results-V4-correlation/Endothelial/k_20/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "C5") %>%  #Y原来是H、C5
  dplyr::select(gs_name, gene_symbol)

#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_POSITIVE_CHEMOTAXIS") %>%
  pull(gene_symbol)
a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "MALAT1") %>%
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 
library(ggpubr)
p <- ggdotchart(a, x = "feature_y", y = "estimate",
                color = "change",                             
                sorting = "ascending",                        
                add = "segments",                             
                add.params = list(color = "lightgray", size = 1.2),
                ggtheme = theme_pubr(),                        
                xlab="",
                dot.size = 2.5)+                              
  labs(title = "GOBP_POSITIVE_CHEMOTAXIS"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of MALAT1") +
  theme(
    axis.line = element_line(size = 0.3, color = "black"),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.text.x = element_text(size = 8,angle = 90,face = "italic",hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black",hjust = 0.5),
    legend.position = "top",
    legend.text=element_text(family="", colour="black", size=8), 
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10))+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.3)+
  scale_color_manual(values = c(low = "#00AFBB",
                                high = "orange"))
ggsave(file.path(output.dir, "bang_MALAT1_chemotaxis.pdf"),
       plot = p,
       height = 2.5,
       width = 4)
#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION") %>%
  pull(gene_symbol)
a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "PCAT19") %>%
  mutate(change = if_else(
    p_value  < 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>% 
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 
library(ggpubr)
p <- ggdotchart(a, x = "feature_y", y = "estimate",
                color = "change",                                
                sorting = "ascending",                        
                add = "segments",                            
                add.params = list(color = "lightgray", size = 1.2),
                ggtheme = theme_pubr(),                        
                xlab="",
                dot.size = 2.5)+                               
  labs(title = "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of PCAT19") +
  theme(
    axis.line = element_line(size = 0.3, color = "black"),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.text.x = element_text(size = 8,angle = 90,face = "italic",hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black",hjust = 0.5),
    legend.position = "top",
    legend.text=element_text(family="", colour="black", size=8), 
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10))+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.3)+
  scale_color_manual(values = c(low = "#00AFBB",
                                high = "orange"))
ggsave(file.path(output.dir, "bang_PCAT19_adhesion.pdf"),
       plot = p,
       height = 2.5,
       width = 4)
#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_VASCULOGENESIS") %>%
  pull(gene_symbol)
a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "MALAT1") %>%
  mutate(change = if_else(
    p_value  < 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 
library(ggpubr)
p <- ggdotchart(a, x = "feature_y", y = "estimate",
                color = "change",                                
                sorting = "ascending",                        
                add = "segments",                             
                add.params = list(color = "lightgray", size = 1.2),
                ggtheme = theme_pubr(),                        
                xlab="",
                dot.size = 2.5)+                             
  labs(title = "GOBP_VASCULOGENESIS"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of MALAT1") +
  theme(
    axis.line = element_line(size = 0.3, color = "black"),
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.text.x = element_text(size = 8,angle = 90,face = "italic",hjust = 1,vjust = 0.5),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black",hjust = 0.5),
    legend.position = "top",
    legend.text=element_text(family="", colour="black", size=8), 
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10))+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.3)+
  scale_color_manual(values = c(low = "#00AFBB",
                                high = "orange"))
p  
ggsave(file.path(output.dir, "bang_MALAT1_vasculogensis.pdf"),
       plot = p,
       height = 2.5,
       width = 4)

##
selected_genes <- c("CYTOR","CD99","GAS6","ANXA1","QKI","CEACAM1","HDAC7","EMP2","VEGFC",
                    "HEY1","HEY2","LEF1","HMGB1",'CX3CL1',"HLA-DRB5","PCAT19",
                    "MALAT1","CD24","CD74","CXCL12","S100A8","SELE","HLA-DRB5",
                    "S100A8","ETS1","HES1")
seurat_obj <- readRDS("~/results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <-
  readRDS("~/results-V4-correlation/Endothelial/k_20/correlation/GO:BP.scores.rds")
scores[1:2, 1:2]
if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]

#
P <- expr %>%
  ggplot(aes(x = MALAT1, y = CXCL12)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MALAT1, y = CXCL12),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  NoLegend()
library(ggExtra)
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "pink"),
                 yparams = list(fill = "#FAE6BE"))
ggsave(file.path(output.dir, "san_MALAT1_CXCL12.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = PCAT19, y = CD74)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = PCAT19, y = CD74),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  NoLegend()
library(ggExtra)
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "pink"),
                 yparams = list(fill = "#FAE6BE"))
ggsave(file.path(output.dir, "san_PCAT19_CD74.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = MALAT1, y = QKI)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MALAT1, y = QKI),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  NoLegend()
library(ggExtra)
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "pink"),
                 yparams = list(fill = "#FAE6BE"))
ggsave(file.path(output.dir, "san_MALAT1_QKI.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = MALAT1, y = HDAC7)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#756bb1", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MALAT1, y = HDAC7),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8, face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, face = "italic", color = "black")) +
  NoLegend()
library(ggExtra)
P2 <- ggMarginal(P, type = "density", 
                 xparams = list(fill = "pink"),
                 yparams = list(fill = "#FAE6BE"))
ggsave(file.path(output.dir, "san_MALAT1_HDAC7.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)

##
adata <-
  readRDS("~/results-V4-correlation/Endothelial/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Endothelial")
# Load package
pacman::p_load(Seurat)
pacman::p_load(tidyverse)

adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

# 
x <- "MALAT1"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-POSITIVE-CHEMOTAXIS",
  group.by = "MALAT1_group",
  assay = "RNA")
p
#绘制boxplot
a <- p$data
colnames(a) <- c("way","ident")
boxplot <- ggpaired(a, 
                    x = "ident", y ="way", color="ident",
                    legend.title=x, 
                    palette = "jco",
                    line.color = "gray", line.size = 0.4,
                    short.panel.labs = FALSE)+
  labs(title = "MALAT1",
       y = "GOBP_POSITIVE_CHEMOTAXIS" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        aspect.ratio = 1.5,
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=0.5),
        plot.title = element_text(face = "italic"),
        axis.title = element_text(size=9), 
        axis.text = element_text(size=9),
        axis.title.x = element_blank()) +
  theme_cat() +
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 0.086,
                        map_signif_level = F, 
                        textsize = 2.5,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) + 
  NoLegend() 
boxplot
ggsave(file.path(output.dir, "deg_MALAT1_chemotaxis.pdf"),
       plot = boxplot,
       height = 2.5,
       width = 2.5)
# 
x <- "PCAT19"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-POSITIVE-REGULATION-OF-LEUKOCYTE-CELL-CELL-ADHESION",
  group.by = "PCAT19_group",
  assay = "RNA")
p
#绘制boxplot
a <- p$data
colnames(a) <- c("way","ident")
boxplot <- ggpaired(a, 
                    x = "ident", y ="way", color="ident",
                    legend.title=x, 
                    palette = "jco",
                    line.color = "gray", line.size = 0.4,
                    short.panel.labs = FALSE)+
  labs(title = "PCAT19",
       y = "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_\nCELL_CELL_ADHESION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        aspect.ratio = 1.5,
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=0.5),
        plot.title = element_text(face = "italic"),
        axis.title = element_text(size=9), 
        axis.text = element_text(size=9),
        axis.title.x = element_blank()) +
  theme_cat() +
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 0.095,
                        map_signif_level = F, 
                        textsize = 2.5,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) + 
  NoLegend() 
boxplot
ggsave(file.path(output.dir, "deg_PCAT19_adhesion.pdf"),
       plot = boxplot,
       height = 2.5,
       width = 2.5)
# 
x <- "MALAT1"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-VASCULOGENESIS",
  group.by = "MALAT1_group",
  assay = "RNA")
p
#绘制boxplot
a <- p$data
colnames(a) <- c("way","ident")

boxplot <- ggpubr::ggpaired(a, 
                            x = "ident", y ="way", color="ident",
                            legend.title=x, 
                            palette = "jco",
                            line.color = "gray", line.size = 0.4,
                            short.panel.labs = FALSE)+
  labs(title = "MALAT1",
       y = "GOBP_VASCULOGENESIS" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        aspect.ratio = 1.5,
        axis.ticks = element_line(colour="black",size=0.4),
        axis.line = element_line(colour="black",size=0.5),
        plot.title = element_text(face = "italic"),
        axis.title = element_text(size=9),
        axis.text = element_text(size=9),
        axis.title.x = element_blank()) +
  theme_cat() +
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 0.14,
                        map_signif_level = F, 
                        textsize = 2.5,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) + 
  NoLegend() 
ggsave(file.path(output.dir, "deg_MALAT1_vasculogenesis.pdf"),
       plot = boxplot,
       height = 2.5,
       width = 2.5)

##
adata <-
  readRDS("~/results-V4-correlation/Endothelial/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Endothelial")
# Load package 
pacman::p_load(Seurat)
pacman::p_load(tidyverse)

adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

# 1、分组MALAT1,图3G_1
x <- "MALAT1"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "CXCL12",  ##LCK  CD3E
  group.by = "MALAT1_group",
  assay = "RNA"
)
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "CXCL12"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "CX3CL1",
  group.by = "MALAT1_group",
  assay = "RNA"
)
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "CX3CL1"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "HMGB1",
  group.by = "MALAT1_group",
  assay = "RNA"
)
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "HMGB1"

f <- rbind(a,b,d)
p<-ggplot(f,
          aes(x=ident, y= value, fill=ident))+
  geom_violin(scale = "area", trim=F, size=0.5,color="white",cex=1,alpha=1)+ 
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  facet_grid(~gene)+
  scale_fill_manual(values=c("orange","#00AFBB"))+ 
  theme_classic(base_size=11, base_family = "",
                base_line_size = 0.3, base_rect_size = 0.3)+
  xlab("")+ylab("MALAT1")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
    axis.ticks = element_line(color="black",size=0.3),
    axis.line = element_line(color="black",size=0.3),
    plot.title = element_text(size = 10, color = "black"))+
  NoLegend()+
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 2.5,
                        map_signif_level = F, 
                        textsize = 3,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) 
p
ggsave(
  file.path(output.dir, "deg_MALAT1_CX3CL1-CXCL12-HMGB1.pdf"),
  plot = p,
  height = 3,
  width = 4.5)

# 
x <- "PCAT19"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "CD74",  ##LCK  CD3E
  group.by = "PCAT19_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "CD74"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "ETS1",
  group.by = "PCAT19_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "ETS1"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "LEF1",
  group.by = "PCAT19_group",
  assay = "RNA")
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "LEF1"

f <- rbind(a,b,d)
p<-ggplot(f,
          aes(x=ident, y= value, fill=ident))+
  geom_violin(scale = "area", trim=F, size=0.5,color="white",cex=1,alpha=1)+ 
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  facet_grid(~gene)+
  scale_fill_manual(values=c("orange","#00AFBB"))+ 
  theme_classic(base_size=11, base_family = "",
                base_line_size = 0.3, base_rect_size = 0.3)+
  xlab("")+ylab("PCAT19")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
    axis.ticks = element_line(color="black",size=0.3),
    axis.line = element_line(color="black",size=0.3),
    plot.title = element_text(size = 10, color = "black"))+
  NoLegend()+
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 4.25,
                        map_signif_level = F, 
                        textsize = 3,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) 
p
ggsave(
  file.path(output.dir, "deg_PCAT19_CD74-ETS1-LEF1.pdf"),
  plot = p,
  height = 3,
  width = 4.5)

# 
x <- "MALAT1"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "EMP2",  ##LCK  CD3E
  group.by = "MALAT1_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "EMP2"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "QKI",
  group.by = "MALAT1_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "QKI"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "HDAC7",
  group.by = "MALAT1_group",
  assay = "RNA")
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "HDAC7"

f <- rbind(a,b,d)
p<-ggplot(f,
          aes(x=ident, y= value, fill=ident))+
  geom_violin(scale = "area", trim=F, size=0.5,color="white",cex=1,alpha=1)+ 
  geom_boxplot(width=0.1,position = position_dodge(0.9),color="white")+
  facet_grid(~gene)+
  scale_fill_manual(values=c("orange","#00AFBB"))+ 
  theme_classic(base_size=11, base_family = "",
                base_line_size = 0.3, base_rect_size = 0.3)+
  xlab("")+ylab("MALAT1")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
    axis.ticks = element_line(color="black",size=0.3),
    axis.line = element_line(color="black",size=0.3),
    plot.title = element_text(size = 10, color = "black"))+
  NoLegend()+
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 1.8,
                        map_signif_level = F, 
                        textsize = 3,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) 
p
ggsave(
  file.path(output.dir, "deg_MALAT1_EMP2-HDAC7-QKI.pdf"),
  plot = p,
  height = 3,
  width = 4.5)

##
PCAT19_adata <- readRDS("~/results-V4-correlation/Endothelial/k_20/deg_PCAT19/C5_gsea.rds")
data <- PCAT19_adata@result
p1 <- PCAT19_adata |>
  cat_gseaplot(
    c("GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION"),
    subplots = c(1, 2, 3),
    color = c("#66C2A5"),#颜色
    title = "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION" |>
      str_replace_all("GOBP_", "") |>
      str_replace_all("_", " ") |>
      str_to_sentence(),
    pvalue_table = T)
p1
ggsave(file.path(output.dir, "gsea_PCAT19_adhesion.pdf"),
       plot = p1,
       height = 2.2,
       width = 3.2)

##
library(ggalluvial)
library(ggplot2)
library(dplyr)
rt <- read.csv("~/results-V13-correlation-multi-lncRNA/End/landscape_Endothelial.csv")
mycol <- c('#ffc1c1',"#FCCDE5","#FFB9B9","#FDB462","#FFBE7A","#f8ac8c","#82B0D2",
           "#66C2A5","#FDE0DF","#CCEBC5","#FDE0DF","#CCEBC5","#8ECFC9","#FDE0DF")
rt$value <- rt$value %>% stringr::str_wrap(width = 20)
p <- ggplot(rt, aes(x = variable, y = Freq,
                    stratum = value, alluvium = flow, fill = value))+ 
  geom_stratum()+
  geom_text(stat='stratum',infer.label=TRUE,size=2.5)+ 
  geom_flow(alpha = .6,aes.flow = "forward")+ 
  scale_fill_manual(values=mycol)+ 
  scale_x_discrete(labels =c('lncRNA','mRNA','pathway'))+
  scale_y_continuous(expand=c(0,0))+ 
  labs(x ='' ,y ='')+
  theme(legend.position = 'none',
        axis.line = element_line(),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())
p
ggsave(file.path(output.dir, "flow.pdf"),
       plot = p,
       height = 3.5,
       width = 4.5)  
