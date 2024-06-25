#### Load packages ----
library(Seurat)
source("~/share/custom_function.R")
source("~/share/custom_plot_function.R")
output.dir <- paste0("~/results-V24-Figure/Fig6") 
selected_cell_type <- "Macrophage"

##
selected_genes <- c("FTX","MALAT1","NEAT1","SNHG8","SNHG9","SNHG29","GAS5","DANCR")
adata <- read_csv("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value < 0.05)
selected_features <- c(
  "GOBP_CELL_CYCLE",
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_REGULATION_OF_MACROPHAGE_PROLIFERATION",
  "GOBP_INTERLEUKIN_12_PRODUCTION",
  "GOBP_INTERLEUKIN_23_PRODUCTION",
  "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
  "GOBP_LIPID_TRANSLOCATION",
  "GOBP_REGULATION_OF_LIPID_TRANSPORT",
  "GOBP_VERY_LOW_DENSITY_LIPOPROTEIN_PARTICLE_ASSEMBLY",
  "GOBP_CHOLESTEROL_EFFLUX",
  "GOBP_REVERSE_CHOLESTEROL_TRANSPORT",
  "GOBP_CHOLESTEROL_CATABOLIC_PROCESS"
)

filtered_adata <- filtered_adata %>% 
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features)

filtered_adata$feature_y <- filtered_adata$feature_y |> 
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

filtered_adata$feature_y <- factor(filtered_adata$feature_y,
                                   levels = selected_features |> 
                                     str_replace_all("GOBP_", "") |>
                                     str_replace_all("_", " ") |>
                                     str_to_sentence())

p <- ggplot(filtered_adata,aes(x = feature_x, y = feature_y, 
                               fill = estimate)) +
  geom_point(aes(size = -log10(p_value)), shape = 21, colour = "black") +
  scale_fill_gradientn(colours = c("#674AA6","white","#CF312C")) + 
  scale_size_continuous(range = c(0,6))+ 
  # coord_fixed() + 
  theme_bw()+ 
  theme(legend.position = "right",
        # title = element_text(size = 8, color = "black",angle = 90), 
        legend.title = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 7, color = "black"),
        panel.spacing.x = unit(0, "pt"),
        panel.grid = element_line(size = 0.2, color = "lightgrey"),
        axis.text.x = element_text(size = 8, color = "black",face = "italic",hjust = 1,vjust = 0.5,angle = 90),#
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, color = "black")) +
  # scale_y_discrete(lables=function(filtered_adata), str_wrap(filtered_adata, width=5))
  scale_y_discrete(labels=function(x) str_wrap(x, width=30))
p 
ggsave(file.path(output.dir, "qipao_pathway_cor.pdf"),
       plot = p,
       height = 3.5,
       width = 4.5)

##
selected_genes <- c("FTX","SNHG29","KCNQ1OT1","OTUD6B-AS1","SNHG8","GAS5")
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
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = GOBP_MACROPHAGE_ACTIVATION,
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_macrophage_activation_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = `SNHG8`, y = GOBP_MACROPHAGE_ACTIVATION)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "GOBP_MACROPHAGE_ACTIVATION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = 0.39, P = 5.1e-18")
ggsave(file.path(output.dir, "san_SNHG8_macrophage_activation.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = GOBP_RESPONSE_TO_INTERFERON_GAMMA,
    #GOCC_IMMUNOGLOBULIN_COMPLEX,
    #GOBP_B_CELL_ACTIVATION,
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_respons_IFNγ_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = `SNHG8`, y = GOBP_RESPONSE_TO_INTERFERON_GAMMA)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), 
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "GOBP_RESPONSE_TO_INTERFERON_GAMMA" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = 0.33, P = 1.3e-13")
p2
ggsave(file.path(output.dir, "san_SNHG8_respons_IFNγ.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `FTX`,
    y = GOBP_CHOLESTEROL_EFFLUX,
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_FTX_cholesterol_efflux_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = `FTX`, y = GOBP_CHOLESTEROL_EFFLUX)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), 
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "CHOLESTEROL_EFFLUX" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = -0.4, P = 6.3e-19")
p2
ggsave(file.path(output.dir, "san_FTX_cholesterol_efflux.pdf"),
       plot = p2,
       height = 2,
       width = 2)

#
p <- expr %>%
  catscatter(
    x = `FTX`,
    y = GOBP_INTERLEUKIN_12_PRODUCTION,
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_FTX_IL12_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = `FTX`, y = GOBP_INTERLEUKIN_12_PRODUCTION)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "GOBP_INTERLEUKIN_12_PRODUCTION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = 0.42, P = 6.3e-21")
p2
ggsave(file.path(output.dir, "san_FTX_IL12.pdf"),
       plot = p2,
       height = 2,
       width = 2)

##
adata <-
  read_csv("~/results-V4-correlation/Macrophage/k_20/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "C5") %>%  #Y原来是H、C5
  dplyr::select(gs_name, gene_symbol)

#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_MACROPHAGE_ACTIVATION") %>%
  pull(gene_symbol)

p <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "SNHG8") %>%
  mutate(change = if_else(
    p_value  < 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%  
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 0.5) +
  geom_point(aes(col = change), size = 0.5) +
  theme_cat() +
  labs(title = "GOBP_MACROPHAGE_ACTIVATION"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of SNHG8") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      face = "italic",
      vjust = 0.5),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(
    low = "#674AA6",
    high = "#CF312C")) +
  NoLegend()
p   
ggsave(file.path(output.dir, "bang_SNHG8_macrophage_activation.pdf"),
       plot = p,
       height = 2.5,
       width = 3)
#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_RESPONSE_TO_INTERFERON_GAMMA") %>%
  pull(gene_symbol)

p <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "SNHG8") %>%
  mutate(change = if_else(
    p_value  < 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  
  top_n(20, wt = abs(estimate)) %>%  
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 0.5) +
  geom_point(aes(col = change), size = 0.5) +
  theme_cat() +
  labs(title = "GOBP_RESPONSE_TO_INTERFERON_GAMMA"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of SNHG8") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      face = "italic",
      vjust = 0.5),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(
    low = "#674AA6",
    high = "#CF312C")) +
  NoLegend()
p   
ggsave(file.path(output.dir, "bang_SNHG8_respons_IFNγ.pdf"),
       plot = p,
       height = 2.5,
       width = 3)
#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_CHOLESTEROL_EFFLUX") %>%
  pull(gene_symbol)

p <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "FTX") %>%
  mutate(change = if_else(
    p_value  < 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  #意思是过滤掉不相关的基因
  top_n(20, wt = abs(estimate)) %>%  #只展示了top40个基因，如果想要看全部的，就看p后面的filtered_adata
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 0.5) +
  geom_point(aes(col = change), size = 0.5) +
  theme_cat() +
  labs(title = "GOBP_CHOLESTEROL_EFFLUX"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of FTX") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      face = "italic",
      vjust = 0.5),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(
    low = "#674AA6",
    high = "#CF312C")) +
  NoLegend()
p   
ggsave(file.path(output.dir, "bang_FTX_cholesterol_efflux.pdf"),
       plot = p,
       height = 2.5,
       width = 3)

#
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_INTERLEUKIN_12_PRODUCTION") %>%
  pull(gene_symbol)

p <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "FTX") %>%
  mutate(change = if_else(
    p_value  < 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation")) %>%
  filter(change != "Uncorrelation") %>%  #意思是过滤掉不相关的基因
  top_n(20, wt = abs(estimate)) %>%  #只展示了top40个基因，如果想要看全部的，就看p后面的filtered_adata
  mutate(feature_y = fct_reorder(feature_y, dplyr::desc(estimate))) %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,
    yend = estimate,
    color = change),
    size = 0.5) +
  geom_point(aes(col = change), size = 0.5) +
  theme_cat() +
  labs(title = "GOBP_INTERLEUKIN_12_PRODUCTION"|>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       y = "Correlation of FTX") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      face = "italic",
      vjust = 0.5),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)) +
  scale_color_manual(values = c(
    low = "#674AA6",
    high = "#CF312C")) +
  NoLegend()
p   
ggsave(file.path(output.dir, "bang_FTX_IL12.pdf"),
       plot = p,
       height = 2.5,
       width = 3)

##
selected_genes <- c("TNF","PTPRC","FOXP1","IL4R",
                    "SNHG8","JUN","LDLR",
                    "FTX","APOE","ABCG1","APOC1","IL10","PLTP",
                    "TNF","IRF1","GBP3","CXCL16")
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
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = JUN,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_JUN_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = SNHG8, y = JUN)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "JUN", 
       title = "R = 0.41, P = 1e-20")
p2
ggsave(file.path(output.dir, "san_SNHG8_JUN.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = TNF,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_TNF_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = SNHG8, y = TNF)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "TNF", 
       title = "R = 0.37, P = 6.1e-17")
p2
ggsave(file.path(output.dir, "san_SNHG8_TNF.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = PTPRC,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_PTPRC_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = SNHG8, y = PTPRC)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "PTPRC", 
       title = "R = 0.48, P = 1.3e-28")
p2
ggsave(file.path(output.dir, "san_SNHG8_PTPRC.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = GBP3,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_GBP3_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = SNHG8, y = GBP3)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "GBP3", 
       title = "R = 0.38, P = 2.8e-17")
p2
ggsave(file.path(output.dir, "san_SNHG8_GBP3.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = IRF1,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_IRF1_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = SNHG8, y = IRF1)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "IRF1", 
       title = "R = 0.35, P = 3.2e-15")
p2
ggsave(file.path(output.dir, "san_SNHG8_IRF1.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `SNHG8`,
    y = CXCL16,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_SNHG8_CXCL16_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = SNHG8, y = CXCL16)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "CXCL16", 
       title = "R = 0.26, P = 7.1e-09")
p2
ggsave(file.path(output.dir, "san_SNHG8_CXCL16.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `FTX`,
    y = ABCG1,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_FTX_ABCG1_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = FTX, y = ABCG1)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "ABCG1", 
       title = "R = -0.4, P = 4.2e-19")
p2
ggsave(file.path(output.dir, "san_FTX_ABCG1.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `FTX`,
    y = APOE,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_FTX_APOE_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = FTX, y = APOE)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "APOE", 
       title = "R = -0.41, P = 1.4e-20")
p2
ggsave(file.path(output.dir, "san_FTX_APOE.pdf"),
       plot = p2,
       height = 2,
       width = 2)
#
p <- expr %>%
  catscatter(
    x = `FTX`,
    y = PLTP,  
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(file.path(output.dir, "san_FTX_PLTP_pre.pdf"),
       plot = p,
       height = 2,
       width = 2)
p2 <- expr %>%
  ggplot(aes(x = FTX, y = PLTP)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#6b76ae",
                        fill = "#e5a323",
                        size = 0.5,
                        alpha = 0.2) +
  theme_cat() +
  theme(#legend.position = "top",
    # title = element_text(size = 8, color = "black",angle = 90), #是legend的
    aspect.ratio = 1, legend.margin = margin(l=-8),
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 7, color = "black"),
    panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 8,face = "italic", color = "black"),
    axis.title.y = element_text(size = 8,face = "italic", color = "black"),
    # axis.text.x = element_text(size = 8, color = "black",face = "italic"),#,angle = 45
    # axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(size = 8, color = "black"))+
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_viridis_c() +
  labs(y = "PLTP", 
       title = "R = -0.43, P = 1.9e-22")
p2
ggsave(file.path(output.dir, "san_FTX_PLTP.pdf"),
       plot = p2,
       height = 2,
       width = 2)
##
adata <-
  readRDS("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Macrophage")
#### Load package
pacman::p_load(Seurat)
pacman::p_load(tidyverse)
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

library(ggunchained)
library(ggsci)
library(Rmisc)
library(ggpubr)
#
x <- "FTX"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("FTX", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-CHOLESTEROL-EFFLUX",
  group.by = "FTX_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$way <- "GOBP_CHOLESTEROL_EFFLUX"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
#画图
Data_summary <- summarySE(a, measurevar="value", groupvars=c("ident","way"))
P1 <- ggplot(data=a, aes(x = way, y = value,fill = ident)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = Data_summary,aes(x = way, y = value),pch=19,
             position=position_dodge(0.5),size=1)+ 
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family="",size=10,colour="black"), 
        axis.title.x = element_text(family="",size = 10,face = "italic"), 
        axis.title.y = element_text(family="",size = 10),
        axis.ticks.x = element_line(colour="black",size=0.3),
        axis.line.x = element_line(colour="black",size=0.3),
        axis.line.y = element_line(colour="black",size=0.3),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(family="", colour="black", size=10), 
        legend.title=element_text(family="", colour="black", size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Cholesterol efflux")+xlab("FTX")+ 
  stat_compare_means(aes(group = ident),
                     label = "p.format",
                     label.y = 0.21,
                     method = "wilcox.test",
                     hide.ns = T)
P1
ggsave(file.path(output.dir, "Vlnplot_FTX_Cholesterol_efflux.pdf"),
       plot = P1,
       height = 2.5,
       width = 3)

#
x <- "FTX"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("FTX", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-INTERLEUKIN-12-PRODUCTION",
  group.by = "FTX_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$way <- "GOBP_INTERLEUKIN_12_PRODUCTION"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
#画图
Data_summary <- summarySE(a, measurevar="value", groupvars=c("ident","way"))
P1 <- ggplot(data=a, aes(x = way, y = value,fill = ident)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = Data_summary,aes(x = way, y = value),pch=19,
             position=position_dodge(0.5),size=1)+ 
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family="",size=10,colour="black"), 
        axis.title.x = element_text(family="",size = 10,face = "italic"), 
        axis.title.y = element_text(family="",size = 10),
        axis.ticks.x = element_line(colour="black",size=0.3),
        axis.line.x = element_line(colour="black",size=0.3),
        axis.line.y = element_line(colour="black",size=0.3),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(family="", colour="black", size=10),
        legend.title=element_text(family="", colour="black", size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Interleukin 12 production")+xlab("FTX")+ 
  stat_compare_means(aes(group = ident),
                     label = "p.format",
                     label.y = 0.13,
                     method = "wilcox.test",
                     hide.ns = T)
P1
ggsave(file.path(output.dir, "Vlnplot_FTX_IL12.pdf"),
       plot = P1,
       height = 2.5,
       width = 3)

#
x <- "SNHG8"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("SNHG8", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)
p2 <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-MACROPHAGE-ACTIVATION",
  group.by = "SNHG8_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$way <- "GOBP_MACROPHAGE_ACTIVATION"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
#画图
Data_summary <- summarySE(b, measurevar="value", groupvars=c("ident","way"))
P1 <- ggplot(data=b, aes(x = way, y = value,fill = ident)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = Data_summary,aes(x = way, y = value),pch=19,
             position=position_dodge(0.5),size=1)+ 
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family="",size=10,colour="black"), 
        axis.title.x = element_text(family="",size = 10,face = "italic"), 
        axis.title.y = element_text(family="",size = 10),
        axis.ticks.x = element_line(colour="black",size=0.3),
        axis.line.x = element_line(colour="black",size=0.3),
        axis.line.y = element_line(colour="black",size=0.3),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(family="", colour="black", size=10),
        legend.title=element_text(family="", colour="black", size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Macrophage activation")+xlab("SNHG8")+ 
  stat_compare_means(aes(group = ident),
                     label = "p.format",
                     label.y = 0.17,
                     method = "wilcox.test",
                     hide.ns = T)
P1
ggsave(file.path(output.dir, "Vlnplot_SNHG8_Macrophage_activation.pdf"),
       plot = P1,
       height = 2.5,
       width = 3)
#
x <- "SNHG8"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("SNHG8", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)
p3 <- VlnPlot(
  sub_seurat_obj,
  features = "GOBP-RESPONSE-TO-INTERFERON-GAMMA",
  group.by = "SNHG8_group",
  assay = "RNA")
p3
c <- p3$data
colnames(c) <- c("value","ident")
c$cell <- rownames(c)
c$way <- "GOBP_RESPONSE_TO_INTERFERON_GAMMA"|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

#画图
Data_summary <- summarySE(c, measurevar="value", groupvars=c("ident","way"))
P1 <- ggplot(data=c, aes(x = way, y = value,fill = ident)) + 
  geom_split_violin(trim=FALSE,color="white") + 
  geom_point(data = Data_summary,aes(x = way, y = value),pch=19,
             position=position_dodge(0.5),size=1)+ 
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(family="",size=10,colour="black"), 
        axis.title.x = element_text(family="",size = 10,face = "italic"), 
        axis.title.y = element_text(family="",size = 10),
        axis.ticks.x = element_line(colour="black",size=0.3),
        axis.line.x = element_line(colour="black",size=0.3),
        axis.line.y = element_line(colour="black",size=0.3),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), 
        legend.text=element_text(family="", colour="black", size=10), 
        legend.title=element_text(family="", colour="black", size=10),
        panel.grid.major = element_blank(),   
        panel.grid.minor = element_blank())+  
  ylab("Response to interferon gama")+xlab("SNHG8")+ 
  stat_compare_means(aes(group = ident),
                     label = "p.format",
                     label.y = 0.2,
                     method = "wilcox.test",
                     hide.ns = T)
P1
ggsave(file.path(output.dir, "Vlnplot_SNHG8_response_gama.pdf"),
       plot = P1,
       height = 2.5,
       width = 3)

##
adata <-
  readRDS("~/results-V4-correlation/Macrophage/k_20/correlation/GO:BP.scores.rds")
seurat_obj <-
  readRDS("~/results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")

sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Macrophage")
#### Load package
pacman::p_load(Seurat)
pacman::p_load(tidyverse)

adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj

#Fig4G-1
x <- "FTX"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("FTX", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "APOE",  ##LCK  CD3E
  group.by = "FTX_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "APOE"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "PLTP",
  group.by = "FTX_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "PLTP"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "ABCG1",
  group.by = "FTX_group",
  assay = "RNA")
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "ABCG1"

f <- rbind(a,b,d)

library(ggpubr)
p <-  ggplot(data = f,mapping =aes(x=ident,y=value,fill=ident))+
  geom_boxplot(size=.3, 
               width=0.4, 
               outlier.size = 0.1)+
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  facet_grid(~gene,scales='free')+
  theme_bw()+ 
  xlab("")+ylab("FTX")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text.x = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank()  
  )+
  geom_signif(
    comparisons = list(c( "high","low")),
    test = "wilcox.test",
    y_position = 5,
    map_signif_level = F, 
    textsize = 2.5,
    size=0.3,
    color = "black",
    tip_length = c(0.02, 0.02)
  )+  NoLegend()
p
ggsave(file.path(output.dir, "boxplot_FTX_ABCG1_PLTP_APOE.pdf"),
       plot = p,
       height = 2.5,
       width = 3.5)

#
x <- "SNHG8"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("SNHG8", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "JUN",  ##LCK  CD3E
  group.by = "SNHG8_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "JUN"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "TNF",
  group.by = "SNHG8_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "TNF"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "PTPRC",
  group.by = "SNHG8_group",
  assay = "RNA")
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "PTPRC"

f <- rbind(a,b,d)

library(ggpubr)
p <-  ggplot(data = f,mapping =aes(x=ident,y=value,fill=ident))+
  geom_boxplot(size=.3, 
               width=0.4, 
               outlier.size = 0.1)+
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  facet_grid(~gene,scales='free')+
  theme_bw()+ 
  xlab("")+ylab("SNHG8")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    # panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text.x = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank()  
  )+
  geom_signif(
    comparisons = list(c("high", "low")),
    test = "wilcox.test",
    y_position = 3.3,
    map_signif_level = F, 
    textsize = 2.5,
    size=0.3,
    color = "black",
    tip_length = c(0.02, 0.02)
  )+NoLegend()
p
ggsave(file.path(output.dir, "boxplot_SNHG8_JUN_PTPRC_TNF.pdf"),
       plot = p,
       height = 2.5,
       width = 3.5)

#
x <- "SNHG8"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c("SNHG8", "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)

p1 <- VlnPlot(
  sub_seurat_obj,
  features = "GBP3",  ##LCK  CD3E
  group.by = "SNHG8_group",
  assay = "RNA")
p1
a <- p1$data 
colnames(a) <- c("value","ident")
a$cell <- rownames(a)
a$gene <- "GBP3"

p2 <- VlnPlot(
  sub_seurat_obj,
  features = "CXCL16",
  group.by = "SNHG8_group",
  assay = "RNA")
p2
b <- p2$data
colnames(b) <- c("value","ident")
b$cell <- rownames(b)
b$gene <- "CXCL16"

p3 <- VlnPlot(
  sub_seurat_obj,
  features = "IRF1",
  group.by = "SNHG8_group",
  assay = "RNA")
p3
d <- p3$data
colnames(d) <- c("value","ident")
d$cell <- rownames(d)
d$gene <- "IRF1"

f <- rbind(a,b,d)

library(ggpubr)
p <-  ggplot(data = f,mapping =aes(x=ident,y=value,fill=ident))+
  geom_boxplot(size=.3, 
               width=0.4, 
               outlier.size = 0.1)+
  scale_fill_manual(values = c("#E99797", "#877FB2"))+ 
  facet_grid(~gene,scales='free')+
  theme_bw()+
  xlab("")+ylab("SNHG8")+
  theme(
    legend.title = element_text(size =10, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    # panel.spacing.x = unit(0, "pt"), #调整各个分面(x轴分面)之间空隙的距离
    #panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.title.x = element_text(size = 10,face = "italic", color = "black"),
    axis.title.y = element_text(size = 10, color = "black",face = "italic"),
    axis.text.x = element_text(size = 10, color = "black"),#,angle = 45,vjust = 0.5
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 10, color = "black"),
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank()  
  )+
  geom_signif(
    comparisons = list(c("high", "low")),
    test = "wilcox.test",
    y_position = 2.3,
    map_signif_level = F, 
    textsize = 2.5,
    size=0.3,
    color = "black",
    tip_length = c(0.02, 0.02)
  )+NoLegend()
p
ggsave(file.path(output.dir, "boxplot_SNHG8_GBP3_IRF1_CXCL16.pdf"),
       plot = p,
       height = 2.5,
       width = 3.5)

##
library(ggalluvial)
library(ggplot2)
library(dplyr)
rt <- read.csv("~/results-V21-Figure4/landscape_macrophage.csv")
mycol <- c("#FAE6BE","#E6F5D0","#f8ac8c","#CCEBC5","#f8ac8c","#82B0D2","#66C2A5",
           "#8ECFC9","#BEB8DC","#FDE0DF","#FCCDE5","#BEB8DC","#BEB8DC","#FFB9B9")

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
