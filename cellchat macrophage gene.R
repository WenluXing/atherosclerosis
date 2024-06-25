#### Load packages ----
library(Seurat)
library(tidyverse)
source("~/share/custom_function.R")
source("~/share/custom_plot_function.R")
output.dir <- paste0("~/results-V24-Figure/Fig7") 
#### Load data ----
adata <-
  readRDS("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")

sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Macrophage")

pacman::p_load(Seurat)
pacman::p_load(tidyverse)

adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

##
x <- "NEAT1"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("NEAT1", "_group")]]))
sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "CXCL8",  
  group.by = "NEAT1_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "CXCL8"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "CXCL3",
  group.by = "NEAT1_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "CXCL3"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "CXCL2",
  group.by = "NEAT1_group",
  assay = "RNA")
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "CXCL2"

f <- rbind(a,b,d)

library(ggpubr)
p <-  ggplot(data = f,mapping =aes(x=ident,y=value,fill=ident))+
  geom_boxplot(size=.3, 
               width=0.6, 
               outlier.size = 0.1)+
  scale_fill_manual(values = c('#E59CC4',"#00AFBB"))+ 
  facet_grid(~gene,scales='free')+
  theme_bw()+ 
  xlab("")+ylab("NEAT1")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    # panel.spacing.x = unit(0, "pt"),
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank()  
  )+
  geom_signif(
    comparisons = list(c( "high","low")),
    test = "wilcox.test",
    y_position = 3.7,
    map_signif_level = F, 
    textsize = 2.5,
    size=0.3,
    color = "black",
    tip_length = c(0.02, 0.02)
  )+  NoLegend()
p
ggsave(file.path(output.dir, "boxplot_NEAT1_CXCL8-CXCL3-CXCL2.pdf"),
       plot = p,
       height = 3,
       width = 3)

##
selected_genes <- c("CXCL2","CXCL3","CXCL8","NEAT1")
seurat_obj <- readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <-
  readRDS("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.scores.rds")
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
  ggplot(aes(x = NEAT1, y = CXCL2)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color="#C388FE",fill="#D1B6E1")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = NEAT1, y = CXCL2),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "CXCL2")+
  NoLegend()
P
ggsave(file.path(output.dir, "Fig7G_san_NEAT1_CXCL2.pdf"),
       plot = P,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = NEAT1, y = CXCL3)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color="#0AA1FF",fill="#a5dff9")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = NEAT1, y = CXCL3),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "CXCL3")+
  NoLegend()
P
ggsave(file.path(output.dir, "san_NEAT1_CXCL3.pdf"),
       plot = P,
       height = 2.5,
       width = 2.5)
#
P <- expr %>%
  ggplot(aes(x = NEAT1, y = CXCL8)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) +
  geom_smooth(method = "lm", formula = y~x, color="#77C034",fill="#C5E99B")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = NEAT1, y = CXCL8),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8,face = "italic", color = "black")) +
  labs(y = "CXCL8")+
  NoLegend()
P
ggsave(file.path(output.dir, "san_NEAT1_CXCL8.pdf"),
       plot = P,
       height = 2.5,
       width = 2.5)

##
selected_genes <- c("NEAT1","CCL2")
adata <- read.csv("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.cor.csv")
seurat_obj <- readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")

table(seurat_obj$cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <- readRDS("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.scores.rds")
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
  ggplot(aes(x = NEAT1, y = GOBP_POSITIVE_REGULATION_OF_CHEMOKINE_PRODUCTION)) +
  geom_point(size = 0.5,color="black",alpha = 0.9) + 
  geom_smooth(method = "lm", formula = y~x, color = "orange", fill = "#cbc9e2")+
  ggpubr::stat_cor(method = 'pearson',
                   aes(x = NEAT1, y = GOBP_POSITIVE_REGULATION_OF_CHEMOKINE_PRODUCTION),
                   label.x.npc = "left",
                   label.y.npc = "top",
                   size = 3,
                   color = "black") +
  theme_cat(aspect.ratio = 1) +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(size = 8,face = "italic", color = "black"),
        axis.title.y = element_text(size = 8, color = "black")) +
  labs(y = "GOBP_POSITIVE_REGULATION_OF_\nCHEMOKINE_PRODUCTION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence())+
  NoLegend()
P
ggsave(file.path(output.dir, "san_NEAT1_chemotaxis.pdf"),
       plot = P,
       height = 2.5,
       width = 2.5)

##
adata <-
  readRDS("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")

sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Macrophage")

pacman::p_load(Seurat)
pacman::p_load(tidyverse)

adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

x <- "NEAT1"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("NEAT1", "_group")]]))
sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

a$gene <- "GOBP_POSITIVE_REGULATION_OF_\nCHEMOKINE_PRODUCTION"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

p <- ggplot(data=a, aes(x = ident, y = way,fill = ident)) + 
  geom_violin(trim = T, position = position_dodge(width = 1),scale = "width",lwd=0.3,alpha=0.8)+
  geom_boxplot(position = position_dodge(1),width=.2,lwd=0.3,alpha=0.5)+
  geom_point(size=0,alpha=0) +
  scale_fill_manual(values = c('#E59CC4',"#00AFBB"))+
  theme_bw()+
  theme(panel.grid = element_line(size = 0.2, color = "lightgrey"),
        axis.title = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.ticks = element_line(color="black",size=0.3))+
  NoLegend()+
  ggsignif::geom_signif(comparisons = list(c("high", "low")),
                        test = "wilcox.test",
                        y_position = 0.15,
                        map_signif_level = F, 
                        textsize = 2.5,
                        size=0.3,
                        color = "black",
                        tip_length = c(0.02, 0.02)) +
  facet_wrap(.~gene,scales = "free",nrow = 1)
p

ggsave(file.path(output.dir, "Voln_NEAT1_chemotaxis.pdf"),
       plot = p,
       height = 2.5,
       width = 2)
