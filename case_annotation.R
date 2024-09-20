library(Seurat)
library(tidyverse)
set.seed(717)
merged <- readRDS("~/data/merged_union_filter.rds")

##
merged <- PercentageFeatureSet(merged,
                               pattern = "^MT-",
                               col.name = "percent.mt")
merged <- PercentageFeatureSet(merged, 
                               pattern = "^RP[SL]", 
                               col.name = "percent_ribo")
merged <- PercentageFeatureSet(merged, 
                               pattern = "^HB[^(P)]", 
                               col.name = "percent_hb")
# 
merged <- CellCycleScoring(merged,
                           g2m.features = cc.genes$g2m.genes,
                           s.features = cc.genes$s.genes)

# 根据可视化结果进行过滤
merged <- subset(merged, subset =
                   nFeature_RNA > 200 &
                   nFeature_RNA < 4000 &
                   percent.mt < 10)

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 3000)
merged <- ScaleData(merged)

##
merged <- RunPCA(merged)
ElbowPlot(merged, ndims = 50)

merged <- merged |>
  harmony::RunHarmony(
    group.by.vars = c("donor"),
    dims.use = 1:20,
    reduction = "pca",
    plot_convergence = TRUE
  )

merged <- merged |>
  RunTSNE(dims = 1:20, reduction = "harmony") |>
  RunUMAP(dims = 1:20, reduction = "harmony") |>
  FindNeighbors(dims = 1:20, reduction = "harmony") |>
  FindClusters(resolution = c(0.8))
saveRDS(merged,"/home/xingwl/practice5-Atherosis/noncoding/anno/anno_reduction_20_0.8.rds")

#
DimPlot(merged,label = T,split.by = "orig.ident",ncol = 3,reduction = "umap")
DimPlot(merged,label = F,group.by = "orig.ident",reduction = "tsne")
DimPlot(merged,
        reduction = "umap",
        label = T)
DimPlot(merged,
        reduction = "tsne", # tsne, umap, pca
        label = T)

##
DefaultAssay(merged) <- "RNA"
features <- c("TPM2","MYL9","ACTA2","TAGLN","MYH11","CNN1","TNS1",
              "DKK3","IGFBP2","BGN","GAS6","TNFRSF11B",           
              "C7","LUM","DCN","PTN","FBLN1","CFD",
              "SERPINE2","OMD","DPT",
              "PDGFRB","CCL21","ID4","FABP4","RGS16","AGT","FGF7","CD36",
              "RERGL","NET1",       
              "ACKR1","CLDN5","VWF","PECAM1","EDN1","GJA5",
              
              "LYZ","MARCO","C1QC","C1QB","C1QA","S100A9","S100A8","CD163","APOC1",             
              "MS4A1","CD79B","CD79A",               
              "CD3D","CD3E","TRAC","IL32","IL7R",
              "KLRD1","PRF1","GNLY","NKG7",               
              "IGHM","MZB1","TXNDC5","EAF2",
              "IGHG1","IGHG2",
              "TPSAB1","CPA3",              
              "S100B","PLP1","GPM6B"
)
DotPlot(merged, features = features, col.min = 0) + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
VlnPlot(merged, features = features, ncol = 3, pt.size = 0)
FeaturePlot(
  merged,
  features = c("CNN1","FN1","TCF21","TNFRSF11B"),
  reduction = "tsne",
  order = T,
  ncol = 2,
  min.cutoff = 0
)

markers <- FindAllMarkers(
  merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

Idents(merged) <- merged$seurat_clusters
top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
write.csv(markers,"~/anno/markers_afteranno_20_0.8.csv", quote = F, row.names = F)

DotPlot(merged, features = unique(top.markers$gene),
        col.min = 0) + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

merged <- ScaleData(object = merged,features = rownames(merged))
DoHeatmap(merged, features = markers$gene,size = 2) + NoLegend()
p <- DoHeatmap(merged, features = top.markers$gene)

#
new.cluster.ids <- c(
  "Macrophage", # 00
  "Endothelial",# 01 
  "Fibroblast", # 02  
  "SMC", # 03
  "T cell", # 04 
  "Fibromyocyte",# 05 
  "Pericyte",# 06 
  "B cell", #07 
  "Fibroblast",#08
  "Pericyte",#09 Pericyte 2
  "Plasma cell",#10Plasma cell 1
  "Macrophage",#11
  "Fibroblast",#12 
  "Macrophage",#13 
  "NK cell",#14 
  "Endothelial",#15 
  "Neuron",#16
  "Endothelial",#17 
  "Fibroblast",#18
  "Mast cell",#19 
  "Fibroblast"#20
)
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged,
                       new.cluster.ids)

merged$cell_type <- Idents(merged)

DimPlot(merged,
        reduction = "tsne",
        label = F,
        pt.size = 0.5,
        cols = c("#E64B35FF","#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF",
                 "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF","#A2005699","#9467BDFF","#42B54099")) 
saveRDS(merged,
        file = "~/anno/harmony_anno_20_0.8.rds")
