#### 加载包 ----
library(Seurat)
source("~/share/custom_function.R")
source("~/share/custom_plot_function.R")
output.dir <- paste0("~/results-V24-Figure/Fig5") 
selected_cell_type <- "SMC"

##
selected_genes <- c("EPB41L4A-AS1","SNHG29","MIR100HG","SNHG9","NEAT1","MALAT1","PHF1","SNHG16")
adata2 <- read_csv("~/results-V4-correlation/SMC/k_20/correlation/H.cor.csv")
adata1 <- read_csv("~/results-V4-correlation/SMC/k_20/correlation/CP:KEGG.cor.csv")
adata <- read_csv("~/results-V4-correlation/SMC/k_20/correlation/GO:BP.cor.csv")

filtered_adata <- adata %>%
  filter(p_value < 0.05)
filtered_adata1 <- adata1 %>%
  filter(p_value < 0.05)
filtered_adata2 <- adata2 %>%
  filter(p_value < 0.05)

selected_features <- c(
  "GOBP_CELL_CYCLE",  
  "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION",
  "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION",
  "GOBP_SMOOTH_MUSCLE_CELL_DIFFERENTIATION",
  "GOBP_BLOOD_VESSEL_REMODELING",
  "GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS",
  "GOBP_PROTEOGLYCAN_BIOSYNTHETIC_PROCESS",
  "GOBP_EXTRACELLULAR_MATRIX_CONSTITUENT_SECRETION",
  "GOBP_LIPID_HOMEOSTASIS",   
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "GOBP_MONOCYTE_CHEMOTACTIC_PROTEIN_1_PRODUCTION",
  "HALLMARK_TGF_BETA_SIGNALING", 
  "GOBP_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY"
)

rt <- filtered_adata %>% 
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features)%>%
  mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  ) 
rt1 <- filtered_adata1 %>% 
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features)%>%
  mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  )
rt2 <- filtered_adata2 %>% 
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features)%>%
  mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  )
rt <- rbind(rt,rt1,rt2)

#相关性r值
rt_r <- rt[,c(1,2,4)]
rt_r <- as.data.frame(rt_r)
rt_r <- spread(rt_r,feature_x,estimate)
rownames(rt_r) <- rt_r$feature_y
rt_r <- rt_r[,-1]
rt_r[is.na(rt_r)] <- 0
#p值
rt_p <- rt[,1:3]
rt_p <- as.data.frame(rt_p)
rt_p <- spread(rt_p,feature_x,p_value)
rownames(rt_p) <- rt_p$feature_y
rt_p <- rt_p[,-1]
rt_p[is.na(rt_p)] <- 0

# a$feature_y <- factor(a$feature_y, levels=selected_features)

min(rt_r)
max(rt_r)
bk1 <- c(seq(-0.8, 0, by = 0.01))
bk2 <- c(seq(0,0.8, by = 0.01))

tolower(rownames(rt_r))
rownames(rt_r) <- selected_features |> 
  str_replace_all("GOBP_", "") |>
  str_replace_all("HALLMARK_", "") |>  
  str_replace_all("_", " ") |>
  str_to_sentence()

library(pheatmap)
p <- pheatmap(rt_r,
              scale = "none",
              cellwidth = 20,
              cellheight = 20,
              number_color="white", 
              number_format="%.2e",
              border="white",
              fontsize_number = 8, 
              display_numbers = matrix(ifelse(rt_p < 0.001, "***",
                                              ifelse(rt_p<0.01,"**",
                                                     ifelse(rt_p<0.05,"*",""))), 
                                       nrow(rt_p)),
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = F,treeheight_col = 20,
              cluster_rows = F,treeheight_row = 20,
              col = c(
                colorRampPalette(colors = c("navy", "#FFFFFF"))(length(bk1)),
                colorRampPalette(colors = c("#FFFFFF", "firebrick3"))(length(bk2))
              ),
              legend_breaks = c(-0.6,0,0.6),
              angle_col = c("90"))
p

ggsave(file.path(output.dir, "heatmap_lncRNA_way.pdf"),
       plot = p,
       height = 7,
       width = 7)

##
#
selected_genes <- c("EPB41L4A-AS1","SNHG29","MIR100HG","SNHG9","NEAT1","MALAT1","PHF1","SNHG16")
seurat_obj <- readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")

table(seurat_obj$cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]

scores <- readRDS("~/results-V4-correlation/SMC/k_20/correlation/GO:BP.scores.rds")
if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]

P <- expr %>%
  ggplot(aes(x = SNHG9, y = GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = SNHG9, y = GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P
library(ggExtra)
P2 <-ggExtra::ggMarginal(P, type = "histogram", 
                         xparams=list(fill = "#3498DB"), 
                         yparams = list(fill="#e23620")) 
P2
ggsave(file.path(output.dir, "san_SNHG9-proliferation.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = SNHG9, y = GOBP_SMOOTH_MUSCLE_CELL_MIGRATION)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = SNHG9, y = GOBP_SMOOTH_MUSCLE_CELL_MIGRATION),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram", 
                         xparams=list(fill = "#3498DB"), 
                         yparams = list(fill="#e23620")) 
P2
ggsave(file.path(output.dir, "san_SNHG9-migration.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)

#
selected_genes <- c("EPB41L4A-AS1","SNHG29","MIR100HG","SNHG9","NEAT1","MALAT1","PHF1","SNHG16")
seurat_obj <- readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")

table(seurat_obj$cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <- readRDS("~/results-V4-correlation/SMC/k_20/correlation/H.scores.rds")
scores[1:2, 1:2]

metacell_obj[["H"]] <- CreateAssayObject(scores)

if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]

P <- expr %>%
  ggplot(aes(x = MIR100HG, y = 	HALLMARK_TNFA_SIGNALING_VIA_NFKB)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MIR100HG, y = HALLMARK_TNFA_SIGNALING_VIA_NFKB),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "HALLMARK_TNFA_SIGNALING_VIA_NFKB" |>
         str_replace_all("HALLMARK_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram",  
                         xparams=list(fill = "#3498DB"), 
                         yparams = list(fill="#e23620")) 
P2
ggsave(file.path(output.dir, "san_MIR100HG_HTNF.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)

#
P <- expr %>%
  ggplot(aes(x = MIR100HG, y = HALLMARK_TGF_BETA_SIGNALING)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MIR100HG, y = HALLMARK_TGF_BETA_SIGNALING),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "HALLMARK_TGF_BETA_SIGNALING" |>
         str_replace_all("HALLMARK_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram",  
                         xparams=list(fill = "#3498DB"), 
                         yparams = list(fill="#e23620")) 
P2
ggsave(file.path(output.dir, "san_MIR100HG_HTGFB.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)

##
adata <-
  read_csv("~/results-V4-correlation/SMC/k_20/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value < 0.05)

#
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "C5") %>%  #Y原来是H、C5
  dplyr::select(gs_name, gene_symbol)
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION") %>%
  pull(gene_symbol)

a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "SNHG9") %>%  
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 

library(ggpubr)
p <- a %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 2) +
  theme_cat() +
  labs(title = "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of SNHG9") +
  theme(
    axis.text.x = element_text(hjust = 1,angle = 60),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(low = "#3498DB",high = "#e23620")) +
  NoLegend()+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.4)
p
ggsave(file.path(output.dir, "bang_SNHG9_proliferation.pdf"),
       plot = p,
       height = 2.5,
       width = 4)
#
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "C5") %>%  #Y原来是H、C5
  dplyr::select(gs_name, gene_symbol)
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION") %>%
  pull(gene_symbol)

a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "SNHG9") %>%  
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 

p <- a %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 2) +
  theme_cat() +
  labs(title = "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of SNHG9") +
  theme(
    axis.text.x = element_text(hjust = 1,angle = 60),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(low = "#3498DB",high = "#e23620")) +
  NoLegend()+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.4)
p
ggsave(file.path(output.dir, "bang_SNHG9_migration.pdf"),
       plot = p,
       height = 2.5,
       width = 4)

#
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "H") %>%  #Y原来是H、C5
  dplyr::select(gs_name, gene_symbol)
selected_features <- gene_sets %>%
  filter(gs_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB") %>%
  pull(gene_symbol)

a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "MIR100HG") %>%  
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 

p <- a %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 2) +
  theme_cat() +
  labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB" |>
         str_replace_all("HALLMARK_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of MIR100HG") +
  theme(
    axis.text.x = element_text(hjust = 1,angle = 60),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(low = "#3498DB",high = "#e23620")) +
  NoLegend()+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.4)
p
ggsave(file.path(output.dir, "bang_MIR100HG_HTNFA.pdf"),
       plot = p,
       height = 2.5,
       width = 4)

#
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "H") %>%  #Y原来是H、C5
  dplyr::select(gs_name, gene_symbol)
selected_features <- gene_sets %>%
  filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>%
  pull(gene_symbol)

a <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "MIR100HG") %>%  
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) 

p <- a %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 2) +
  theme_cat() +
  labs(title = "HALLMARK_TGF_BETA_SIGNALING" |>
         str_replace_all("HALLMARK_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of MIR100HG") +
  theme(
    axis.text.x = element_text(hjust = 1,angle = 60),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(low = "#3498DB",high = "#e23620")) +
  NoLegend()+
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.4)
p
ggsave(file.path(output.dir, "bang_MIR100HG_HTGF.pdf"),
       plot = p,
       height = 2.5,
       width = 4)

##
selected_genes <- c("MIR100HG","SNHG9","NEAT1","FOXP1","ADAMTS1","TGFBR2","PDGFRB","PDGFD","MYC",
                    "EGR1","CCL2","TGIF1","CEBPD","KLF4","JUNB","BMPR2","ID1")
seurat_obj <- readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)

#
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]

scores <-
  readRDS("~/results-V4-correlation/SMC/k_20/correlation/GO:BP.scores.rds")
scores[1:2, 1:2]

if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]
P <- expr %>%
  ggplot(aes(x = SNHG9, y = FOXP1)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x =  SNHG9, y = FOXP1),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "FOXP1")+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram",  
                         xparams=list(fill = "#6b76ae"),  
                         yparams = list(fill="#e5a323")) 
P2
ggsave(file.path(output.dir, "san_SNHG9_FOXP1.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = SNHG9, y = PDGFRB)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = SNHG9, y = PDGFRB),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "PDGFRB")+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram",  
                         xparams=list(fill = "#6b76ae"),  
                         yparams = list(fill="#e5a323")) 
P2
ggsave(file.path(output.dir, "san_SNHG9_PDGFRB.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <-
  readRDS("~/results-V4-correlation/SMC/k_20/correlation/H.scores.rds")
scores[1:2, 1:2]
metacell_obj[["H"]] <- CreateAssayObject(scores)

if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]
P <- expr %>%
  ggplot(aes(x = MIR100HG, y = CEBPD)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MIR100HG, y = CEBPD),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "CEBPD")+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram",  
                         xparams=list(fill = "#6b76ae"),  
                         yparams = list(fill="#e5a323")) 
P2
ggsave(file.path(output.dir, "san_MIR100HG_CEBPD.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = MIR100HG, y = TGIF1)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color = "#E64B35CC", fill = "#F39B7FCC")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = MIR100HG, y = TGIF1),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "TGIF1")+
  NoLegend()
P
P2 <-ggExtra::ggMarginal(P, type = "histogram",  
                         xparams=list(fill = "#6b76ae"),  
                         yparams = list(fill="#e5a323")) 
P2
ggsave(file.path(output.dir, "san_MIR100HG_TGIF1.pdf"),
       plot = P2,
       height = 2.5,
       width = 2.5)

##
adata <-
  readRDS("~/results-V4-correlation/SMC/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "SMC")
#### Load package 
pacman::p_load(Seurat)
pacman::p_load(tidyverse)
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj
# FigS5B-1
x <- "SNHG9"
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
  features = "GOBP-SMOOTH-MUSCLE-CELL-PROLIFERATION",
  group.by = "SNHG9_group",
  assay = "RNA")
p1 
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$way <- "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

p <- ggviolin(a, x="ident", y="value", fill = "ident", 
              add = "boxplot", add.params = list(fill="white"))+ 
  scale_color_manual(values = c("#008B4533","#ddc5b2"))+
  scale_fill_manual(values = c("#008B4533","#ddc5b2"))+
  geom_point(aes(col=ident),shape=1,size=1.5,alpha=0.7,
             position = position_jitterdodge(jitter.width = 0.4,
                                             jitter.height = 0,
                                             dodge.width = 1)) +
  stat_summary(fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "white",position =position_dodge(1)) +
  geom_signif(comparisons = list(c( "high","low")),
              test = "wilcox.test",
              y_position = 0.12,#调整显示p值横线的位置
              map_signif_level = F, 
              textsize = 2.5,#文本的大小
              size=0.3,#调整横线的粗细
              color = "black",#调整横线及p值文本的颜色
              tip_length = c(0.03, 0.03))+#调整左右竖线的长短
  theme(strip.text.x = element_text(size=9,color="black"),
        axis.ticks = element_line(colour="black",size=0.3),
        axis.line = element_line(colour="black",size=0.3),
        plot.title = element_text(size = 10, colour="black", hjust = 0.5),
        axis.text = element_text(size=9),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black",face = "italic"))+
  labs(title="GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),x = "",y="SNHG9")+
  NoLegend()+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
p
ggsave(file.path(output.dir, "Vlnplot_SNHG9_proliferation.pdf"),
       plot = p,
       height = 2.5,
       width = 2.6)
#
p2 <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-SMOOTH-MUSCLE-CELL-MIGRATION",
  group.by = "SNHG9_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$way <- "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

p <- ggviolin(b, x="ident", y="value", fill = "ident", 
              add = "boxplot", add.params = list(fill="white"))+ 
  scale_color_manual(values = c("#008B4533","#ddc5b2"))+
  scale_fill_manual(values = c("#008B4533","#ddc5b2"))+
  geom_point(aes(col=ident),shape=1,size=1.5,alpha=0.7,
             position = position_jitterdodge(jitter.width = 0.4,
                                             jitter.height = 0,
                                             dodge.width = 1)) +
  stat_summary(fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "white",position =position_dodge(1)) +
  geom_signif(comparisons = list(c( "high","low")),
              test = "wilcox.test",
              y_position = 0.125,#调整显示p值横线的位置
              map_signif_level = F, 
              textsize = 2.5,#文本的大小
              size=0.3,#调整横线的粗细
              color = "black",#调整横线及p值文本的颜色
              tip_length = c(0.03, 0.03))+#调整左右竖线的长短
  theme(strip.text.x = element_text(size=9,color="black"),
        axis.ticks = element_line(colour="black",size=0.3),
        axis.line = element_line(colour="black",size=0.3),
        plot.title = element_text(size = 10, colour="black", hjust = 0.5),
        axis.text = element_text(size=9),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black",face = "italic"))+
  labs(title="GOBP_SMOOTH_MUSCLE_CELL_MIGRATION"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),x = "",y="SNHG9")+
  NoLegend()+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
p
ggsave(file.path(output.dir, "Vlnplot_SNHG9_migration.pdf"),
       plot = p,
       height = 2.5,
       width = 2.5)
#
adata <-
  readRDS("~/results-V4-correlation/SMC/k_20/correlation/H.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "SMC")
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

x <- "MIR100HG"
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

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "HALLMARK-TNFA-SIGNALING-VIA-NFKB",
  group.by = "MIR100HG_group",
  assay = "RNA"
)
p3
f <- p3$data
colnames(f) <- c("value","ident")
f$cell <- rownames(f)
f$way <- "HALLMARK_TNFA_SIGNALING_VIA_NFKB"|>
  str_replace_all("HALLMARK_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

p <- ggviolin(f, x="ident", y="value", fill = "ident", 
              add = "boxplot", add.params = list(fill="white"))+ 
  scale_color_manual(values = c("#008B4533","#ddc5b2"))+
  scale_fill_manual(values = c("#008B4533","#ddc5b2"))+
  geom_point(aes(col=ident),shape=1,size=1.5,alpha=0.7,
             position = position_jitterdodge(jitter.width = 0.4,
                                             jitter.height = 0,
                                             dodge.width = 1)) +
  stat_summary(fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "white",position =position_dodge(1)) +
  geom_signif(comparisons = list(c( "high","low")),
              test = "wilcox.test",
              y_position = 0.143,#调整显示p值横线的位置
              map_signif_level = F, 
              textsize = 2.5,#文本的大小
              size=0.3,#调整横线的粗细
              color = "black",#调整横线及p值文本的颜色
              tip_length = c(0.03, 0.03))+#调整左右竖线的长短
  theme(strip.text.x = element_text(size=9,color="black"),
        axis.ticks = element_line(colour="black",size=0.3),
        axis.line = element_line(colour="black",size=0.3),
        plot.title = element_text(size = 10, colour="black", hjust = 0.5),
        axis.text = element_text(size=9),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black",face = "italic"))+
  labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB"|>
         str_replace_all("HALLMARK_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),x = "",y="MIR100HG")+
  NoLegend()+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
p
ggsave(file.path(output.dir, "Vlnplot_MIR100HG_HTNFA.pdf"),
       plot = p,
       height = 2.5,
       width = 2.5)
#
p4 <- VlnPlot(
  sub_seurat_obj,
  features = "HALLMARK-TGF-BETA-SIGNALING",
  group.by = "MIR100HG_group",
  assay = "RNA")
p4
d <- p4$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$way <- "HALLMARK_TGF_BETA_SIGNALING"|>
  str_replace_all("HALLMARK_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

p <- ggviolin(d, x="ident", y="value", fill = "ident", 
              add = "boxplot", add.params = list(fill="white"))+ 
  scale_color_manual(values = c("#008B4533","#ddc5b2"))+
  scale_fill_manual(values = c("#008B4533","#ddc5b2"))+
  geom_point(aes(col=ident),shape=1,size=1.5,alpha=0.7,
             position = position_jitterdodge(jitter.width = 0.4,
                                             jitter.height = 0,
                                             dodge.width = 1)) +
  stat_summary(fun.y = "median",geom = "point",shape = 16, 
               size = 2, color = "white",position =position_dodge(1)) +
  geom_signif(comparisons = list(c( "high","low")),
              test = "wilcox.test",
              y_position = 0.143,#调整显示p值横线的位置
              map_signif_level = F, 
              textsize = 2.5,#文本的大小
              size=0.3,#调整横线的粗细
              color = "black",#调整横线及p值文本的颜色
              tip_length = c(0.03, 0.03))+#调整左右竖线的长短
  theme(strip.text.x = element_text(size=9,color="black"),
        axis.ticks = element_line(colour="black",size=0.3),
        axis.line = element_line(colour="black",size=0.3),
        plot.title = element_text(size = 10, colour="black", hjust = 0.5),
        axis.text = element_text(size=9),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black",face = "italic"))+
  labs(title="HALLMARK_TGF_BETA_SIGNALING"|>
         str_replace_all("HALLMARK_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),x = "",y="MIR100HG")+
  NoLegend()+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
p
ggsave(file.path(output.dir, "Vlnplot_MIR100HG_HTGFB.pdf"),
       plot = p,
       height = 2.5,
       width = 2.5)

##
adata <-
  readRDS("~/results-V4-correlation/SMC/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "SMC")
#### Load package 
pacman::p_load(Seurat)
pacman::p_load(tidyverse)
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

# Fig5D-1/2
x <- "SNHG9"
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
p1 <- VlnPlot(sub_seurat_obj,
              features = c("CTNNB1"),
              group.by = "SNHG9_group",
              assay = "RNA")
p2 <- VlnPlot(sub_seurat_obj,
              features = c("FOXP1"),
              group.by = "SNHG9_group",
              assay = "RNA")
p3 <- VlnPlot(sub_seurat_obj,
              features = c("ADAMTS1"),
              group.by = "SNHG9_group",
              assay = "RNA")
p4 <- VlnPlot(sub_seurat_obj,
              features = c("PDGFRB"),
              group.by = "SNHG9_group",
              assay = "RNA")
p5 <- VlnPlot(sub_seurat_obj,
              features = c("PDGFD"),
              group.by = "SNHG9_group",
              assay = "RNA")
p6 <- VlnPlot(sub_seurat_obj,
              features = c("MYC"),
              group.by = "SNHG9_group",
              assay = "RNA")
a <- p1$data
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "CTNNB1"
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "FOXP1"
c <- p3$data
colnames(c) <- c("value","ident")
c$cell <- rownames(c)
c$gene <- "ADAMTS1"
d <- p4$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "PDGFRB"
f <- p5$data
colnames(f) <- c("value","ident")
f$cell <- rownames(f)
f$gene <- "PDGFD"
g <- p6$data
colnames(g) <- c("value","ident")
g$cell <- rownames(g)
g$gene <- "MYC"
h <- rbind(a,b,c,d,f,g)

p1 <- ggplot(h,aes(x=ident,y=value,fill=ident))+
  scale_fill_manual(values = c('#eb4b3a',"#48bad0"))+
  guides(fill=guide_legend(title=x))+
  labs(x = "", y = "SNHG9")+
  geom_boxplot()+ facet_wrap(~gene,nrow =1)+ theme_bw()+
  theme(axis.line = element_line(color = "black",size=0.08), #将x=0轴和y=0轴加粗显示(size=1)
        axis.ticks = element_line(color = "black",size=0.2),
        axis.title.x = element_text(size = 10,face = "italic", color = "black"),
        axis.title.y = element_text(size = 10, color = "black",face = "italic"),
        axis.text.x = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),  #不显示网格线
        panel.grid.minor = element_blank())+
  geom_signif(comparisons = list(c( "high","low")),
              test = "wilcox.test",
              y_position = 2.35,#调整显示p值横线的位置
              map_signif_level = F,
              textsize = 2.3,#文本的大小
              size=0.3,#调整横线的粗细
              color = "black",#调整横线及p值文本的颜色
              tip_length = c(0.03, 0.03))+#调整左右竖线的长短
  NoLegend()
p1
ggsave(file.path(output.dir, "boxplot_SNHG9_ADAMTS1-FOXP1-CTNNB1_PDGFD-PDGFRB-MYC.pdf"),
       plot = p1,
       height = 3,
       width = 6)

#
adata <-
  readRDS("~/results-V4-correlation/SMC/k_20/correlation/H.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "SMC")
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

x <- "MIR100HG"
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
p1 <- VlnPlot(sub_seurat_obj,
              features = c("EGR1"),
              group.by = "MIR100HG_group",
              assay = "RNA")
p2 <- VlnPlot(sub_seurat_obj,
              features = c("CEBPD"),
              group.by = "MIR100HG_group",
              assay = "RNA")
p3 <- VlnPlot(sub_seurat_obj,
              features = c("JUNB"),
              group.by = "MIR100HG_group",
              assay = "RNA")
a <- p1$data
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "EGR1"
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "CEBPD"
c <- p3$data
colnames(c) <- c("value","ident")
c$cell <- rownames(c)
c$gene <- "JUNB"
d <- rbind(a,b,c)

p1 <- ggplot(d,aes(x=ident,y=value,fill=ident))+
  scale_fill_manual(values = c('#eb4b3a',"#48bad0"))+
  guides(fill=guide_legend(title=x))+
  labs(x = "", y = "MIR100HG")+
  geom_boxplot()+ facet_wrap(~gene,nrow =1)+ theme_bw()+
  theme(axis.line = element_line(color = "black",size=0.08), #将x=0轴和y=0轴加粗显示(size=1)
        axis.ticks = element_line(color = "black",size=0.2),
        axis.title.x = element_text(size = 10,face = "italic", color = "black"),
        axis.title.y = element_text(size = 10, color = "black",face = "italic"),
        axis.text.x = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
        axis.text.y = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),  #不显示网格线
        panel.grid.minor = element_blank())+
  geom_signif(comparisons = list(c( "high","low")),
              test = "wilcox.test",
              y_position = 3.5,#调整显示p值横线的位置
              map_signif_level = F,
              textsize = 2.3,#文本的大小
              size=0.3,#调整横线的粗细
              color = "black",#调整横线及p值文本的颜色
              tip_length = c(0.03, 0.03))+#调整左右竖线的长短
  NoLegend()
p1
ggsave(file.path(output.dir, "boxplot_MIR100HG_CEBPD-JUNB-EGR1.pdf"),
       plot = p1,
       height = 3,
       width = 3)

##
#
SNHG9_adata <- readRDS("~/results-V4-correlation/SMC/k_20/SNHG9_deg/C5_gsea.rds")
data <- SNHG9_adata@result
p1 <- SNHG9_adata |>
  cat_gseaplot(
    c("GOBP_SMOOTH_MUSCLE_CELL_MIGRATION"),
    subplots = c(1, 2, 3),
    color = c("#008B4599"),#颜色"#00A087"
    title = "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION" |>
      str_replace_all("GOBP_", "") |>
      str_replace_all("_", " ") |>
      str_to_sentence(),
    pvalue_table = T)
ggsave(file.path(output.dir, "gsea_SNHG9_SMC_Migration.pdf"),
       plot = p1,
       height = 2.4,
       width = 3)
#
p2 <- SNHG9_adata |>
  cat_gseaplot(
    c("GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION"),
    subplots = c(1, 2, 3),
    color = c("#EE0000FF"),#颜色"#00A087"
    title = "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION" |>
      str_replace_all("GOBP_", "") |>
      str_replace_all("_", " ") |>
      str_to_sentence(),
    pvalue_table = T)
ggsave( file.path(output.dir, "gsea_SNHG9_SMC_Proliferation.pdf"),
        plot = p2,
        height = 2.4,
        width = 3)
#
MIR100HG_adata <- readRDS("~/results-V4-correlation/SMC/k_20/MIR100HG_deg/h_gsea.rds")
p3 <- MIR100HG_adata |>
  cat_gseaplot(
    c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    subplots = c(1, 2, 3),
    color = c("#631879FF"),#颜色
    title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB" |>
      str_replace_all("HALLMARK_", "") |>
      str_replace_all("_", " ") |>
      str_to_sentence(),
    pvalue_table = T)
p3
ggsave(file.path(output.dir, "gsea_MIR100HG_TNFA.pdf"),
       plot = p3,
       height = 2.4,
       width = 3)
#
p4 <- MIR100HG_adata |>
  cat_gseaplot(
    c("HALLMARK_TGF_BETA_SIGNALING"),
    subplots = c(1, 2, 3),
    color = c("#FF6600"),#颜色
    title = "HALLMARK_TGF_BETA_SIGNALING" |>
      str_replace_all("HALLMARK_", "") |>
      str_replace_all("_", " ") |>
      str_to_sentence(),
    pvalue_table = T)
p4
ggsave( file.path(output.dir, "gsea_MIR100HG_TGFB.pdf"),
        plot = p4,
        height = 2.4,
        width = 3)
##
library(ggalluvial)
library(ggplot2)
library(dplyr)
rt <- read.csv("~/results-V21-Figure5/landscape_SMC.csv")
mycol <- c("#FFB9B9","#FCCDE5","#FDE0DF","#FDB462","#FFBE7A", "#f8ac8c",
           "#B3DE69","#CCEBC5","#FCCDE5","#E6F5D0","#BEB8DC","#80B1D3",
           "#8ECFC9","#8ECFC9","#8ECFC9","#FCCDE5","#66C2A5","#FCCDE5")

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
ggsave(file.path(output.dir, "=flow.pdf"),
       plot = p,
       height = 4,
       width = 5.5)  

##
library(WGCNA)
library(dplyr)
source("~/share/20220802_scrna-m6A/custom_function.R")
source("~/share/20220802_scrna-m6A/custom_plot_function.R")
output.dir <- paste0("~/results-V24-Figure/Fig5/wgcna_pre2") 

##
##
seurat.obj <- readRDS("~/results-V15-metacell-separate-celltype-gd/NEAT1_gd/SMC_NEAT1_gd_allmarker-813.rds")
expr <- as.data.frame(seurat.obj[["RNA"]]@data)
dim(expr)
expr <- t(expr)

##
seurat.obj@meta.data$cell <- rownames(seurat.obj@meta.data)
phen <- seurat.obj@meta.data[,7:8]
write.csv(phen,file.path(output.dir, "phen.csv"))

phen$NEAT1_high <- phen$NEAT1_group
phen$NEAT1_low <- phen$NEAT1_group
phen$NEAT1_high[phen$NEAT1_high == "high"] <- "1" 
phen$NEAT1_high[phen$NEAT1_high == "low"] <- "0" 
phen$NEAT1_low[phen$NEAT1_low == "high"] <- "0" 
phen$NEAT1_low[phen$NEAT1_low == "low"] <- "1" 
phen <- as.data.frame(phen[,3:4])
phen$NEAT1_high <- as.numeric(phen$NEAT1_high)
phen$NEAT1_low <- as.numeric(phen$NEAT1_low)
write.csv(phen,file.path(output.dir, "phen-V1.csv"))

##
traitColors = numbers2colors(phen, signed = FALSE)
sampleTree2 = hclust(dist(expr), method = "average")
pdf(file.path(output.dir,"Cell dendrogram and phen heatmap.pdf"),width=6,height=4)
P1 <- plotDendroAndColors(sampleTree2,
                          traitColors,
                          groupLabels = names(phen),
                          main = "Cell dendrogram and phen heatmap")
print(P1)
dev.off()

saveRDS(sampleTree2,file.path(output.dir, "sampleTree2.rds"))

##
powers <- c(1:20)
sft <- pickSoftThreshold(expr,
                         powerVector = powers, verbose = 5)
saveRDS(sft,file.path(output.dir,"sft.rds"))

par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",
     main = paste("Scale independence"))

text(sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = cex1,
     col = "red") +
  abline(h = 0.8, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1,
     col = "red") + abline(h = 100, col = "red")

# 
sft$powerEstimate
softPower <- 3
adjacency <- adjacency(expr, power = softPower)
TOM <- TOMsimilarity(adjacency)
saveRDS(TOM,file.path(output.dir,"TOM.rds"))
dissTOM = 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

# 
sizeGrWindow(12, 9)
plot(geneTree,
     xlab = "",
     sub = "",
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)

# 
minModuleSize <- 30

# 
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 0,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
# 
table(dynamicMods)
# 
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# 
sizeGrWindow(8, 6)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList <- moduleEigengenes(expr, colors = dynamicColors) 
saveRDS(MEList,file.path(output.dir,"MEList.rds"))

MEs <- MEList$eigengenes
MEs$MEgrey <- "0" 
MEs$MEblue <- as.numeric(MEs$MEblue)
MEs$MEbrown <- as.numeric(MEs$MEbrown)
MEs$MEgreen <- as.numeric(MEs$MEgreen)
MEs$MEturquoise <- as.numeric(MEs$MEturquoise)
MEs$MEyellow <- as.numeric(MEs$MEyellow)
MEs$MEgrey <- as.numeric(MEs$MEgrey)
# Calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs)
MEDiss <- as.data.frame(MEDiss)
write.csv(MEDiss,file.path(output.dir,"MEDiss.csv"))
MEDiss$MEgrey <- "0" 
MEDiss[6,] <- "0" 

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
pdf(file.path(output.dir,"Clustering of module eigengenes.pdf"), width = 5, height = 4)
P0 <- plot(METree,
           main = "Clustering of module eigengenes",
           xlab = "",
           sub = "")
print(P0)
dev.off()

MEDissThres <- 0.1
abline(h = MEDissThres, col = "red")
merge = mergeCloseModules(expr,
                          dynamicColors,
                          cutHeight = MEDissThres,
                          verbose = 3) #有一丢费时间
saveRDS(merge,file.path(output.dir,"merge.rds"))

mergedColors <- merge$colors

length(table(mergedColors))

mergedMEs <- merge$newMEs

##
sizeGrWindow(12, 9)
pdf(file.path(output.dir,"Dynamic Tree Cut_Merged dynamic.pdf"), width = 6, height = 4)
P2 <- plotDendroAndColors(geneTree,
                          cbind(dynamicColors, mergedColors),
                          c("Dynamic Tree Cut", "Merged dynamic"),
                          dendroLabels = FALSE,
                          hang = 0.03,
                          addGuide = TRUE,
                          guideHang = 0.05)
print(P2)
dev.off()

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))

moduleLabels <- match(moduleColors, colorOrder) - 1

MEs <- mergedMEs

length(table(moduleColors))
table(moduleColors)

##
nGenes <- ncol(expr)
nSamples <- nrow(expr)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(expr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <-
  cor(MEs, phen, use = "p")
#
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 
                                      nSamples)

# 
sizeGrWindow(15, 20)
textMatrix =  paste(signif(moduleTraitCor, 2),
                    "\n(",
                    signif(moduleTraitPvalue, 1),
                    ")",
                    sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 8, 2, 2))

# Display the correlation values within a heatmap plot
pdf(file.path(output.dir,"Module-trait relationships.pdf"), width = 3.5, height = 5)
P3 <- labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(phen),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors.lab.x = "black",
  colors.lab.y = "black",
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  cex.lab = 0.6,
  zlim = c(-0.8, 0.8),
  main = paste("Module-trait relationships"))
print(P3)
dev.off()

# 
months <- as.data.frame(phen$NEAT1_high)

names(months) = "NEAT1_high"
# names (colors) of the_ modules
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(expr, 
                                          MEs, 
                                          use = "p"))

MMPvalue <-
  as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(expr, months, use = "p"))

GSPvalue <-
  as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <-
  paste("GS.", names(months), sep = "")

names(GSPvalue) <- paste("p.GS.", names(months), sep = "")

# 
module = "yellow"

column = match(module, modNames)

moduleGenes = moduleColors == module

##
sizeGrWindow(7, 7)
par(mfrow = c(1, 1))
pdf(file.path(output.dir,"Module membership vs. gene significance.pdf"), width = 5, height = 5)
P4 <- verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for proerythroblast",
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = module)
print(P4)
dev.off()

##
moduleColors <- mergedColors
module = "yellow"
# Select module probes
probes <- colnames(expr) 
inModule <- (moduleColors == module)
table(inModule)
modProbes <- probes[inModule]

head(modProbes)
length(modProbes)
modGenes <- modProbes

##
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

go <- enrichGO(
  gene = modGenes,
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1)
a <- go@result
write.csv(a,file.path(output.dir,"GO_BP_yellow.csv"))

pathways <- c("myofibril assembly",
              "actin-myosin filament sliding",
              "actomyosin structure organization",
              "muscle system process",
              "muscle cell differentiation",
              "muscle contraction",
              "actin filament-based movement",
              "actin filament organization")
a <- a[a$Description %in% pathways,]

a$Description <- a$Description |>
  str_to_sentence()
##
p <- a %>%
  mutate(Description = fct_reorder(Description,-log10(pvalue))) %>% 
  ggplot(aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", fill = "#FE8D3C") +
  geom_vline(xintercept = 2, linetype = 2) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "-log10(P value)", title = "MEyellow") +
  theme_cat() +
  theme(aspect.ratio = 1, axis.title.y = element_blank(),
        axis.title = element_text(size = 8, color = "black"),
        axis.line = element_line(colour="black",size=0.1),
        axis.ticks = element_line(colour="black",size=0.1))
p
ggsave(file.path(output.dir, "GO_BP_yellow.pdf"),
       plot = p,
       height = 2.5,
       width = 3.5)

##
MEs = moduleEigengenes(expr, moduleColors)$eigengenes
weight = as.data.frame(phen[,c("NEAT1_high","NEAT1_low")])
MET = orderMEs(cbind(MEs, weight))
sizeGrwindow(5,7.5)
par(cex = 0.9)
pdf(file.path(output.dir,"Eigengene dendrogram_adjacency heatmap.pdf"), width = 4, height = 5)
P6 <- plotEigengeneNetworks(MET,"", marDendro =c(0,4,1,2),marHeatmap =c(3,4,1,2),
                            cex.lab = 0.8, xLablesAngle=90)
print(P6)
dev.off()

sizeGrwindow(6,6)
par(cex = 1.0)
pdf(file.path(output.dir,"Eigengene dendrogram.pdf"), width = 4, height = 2.5)
P7 <- plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                            plotHeatmaps = FALSE)
print(P7)
dev.off()

par(cex = 1.0)
pdf(file.path(output.dir,"Eigengene adjacency heatmap.pdf"), width = 4, height = 3)
P8 <- plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                            plotDendrograms=FALSE, xLabelsAngle=90)
print(P8)
dev.off()
