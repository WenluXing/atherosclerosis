##
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(hdf5r)
library(Seurat)

##
datadir <- "~/20240103_Atherosis/CONTROL/data/ATVB2021_use/cardiac_arteries_processed_data"
P1_CA <- Read10X(data.dir = paste0(datadir,"/patient1_CA/"))
P1_CA <- CreateSeuratObject(counts = P1_CA, project = "patient1_CA", min.cells = 3, min.features = 200)
P1_CA

head(P1_CA)
P1_CA$donor <- "patient1"
P1_CA$region <- "Coronary"
P1_CA$diagnosis <- "end-stage heart failure but no discernible atherosclerotic lesions"

P2_CA1 <- Read10X(data.dir = paste0(datadir,"/patient2_CA1/"))
P2_CA1 <- CreateSeuratObject(counts = P2_CA1, project = "patient2_CA1", min.cells = 3, min.features = 200)
P2_CA1

head(P2_CA1)
P2_CA1$donor <- "patient2"
P2_CA1$region <- "Coronary"
P2_CA1$diagnosis <- "End-stage heart failure but no discernible atherosclerotic lesions"

P2_CA2 <- Read10X(data.dir = paste0(datadir,"/patient2_CA2/"))
P2_CA2 <- CreateSeuratObject(counts = P2_CA2, project = "patient2_CA2", min.cells = 3, min.features = 200)
P2_CA2

head(P2_CA2)
P2_CA2$donor <- "patient2"
P2_CA2$region <- "Coronary"
P2_CA2$diagnosis <- "End-stage heart failure but no discernible atherosclerotic lesions"

P3_CA1 <- Read10X(data.dir = paste0(datadir,"/patient3_CA1/"))
P3_CA1 <- CreateSeuratObject(counts = P3_CA1, project = "patient3_CA1", min.cells = 3, min.features = 200)
P3_CA1

head(P3_CA1)
P3_CA1$donor <- "patient3"
P3_CA1$region <- "Coronary"
P3_CA1$diagnosis <- "End-stage heart failure but no discernible atherosclerotic lesions"

P3_CA2 <- Read10X(data.dir = paste0(datadir,"/patient3_CA2/"))
P3_CA2 <- CreateSeuratObject(counts = P3_CA2, project = "patient3_CA2", min.cells = 3, min.features = 200)
P3_CA2

head(P3_CA2)
P3_CA2$donor <- "patient3"
P3_CA2$region <- "Coronary"
P3_CA2$diagnosis <- "End-stage heart failure but no discernible atherosclerotic lesions"

merged <- merge(P1_CA,c(P2_CA1,P2_CA2,P3_CA1,P3_CA2))
table(merged$orig.ident)
table(merged$donor)
head(merged)
saveRDS(merged,"1_merged_unqc.rds")


##
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 20)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
merged

merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)；
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
ElbowPlot(merged, ndims = 50)
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.3)

merged <- RunUMAP(merged, dims = 1:30)
merged <- RunTSNE(merged, dims = 1:30)

DimPlot(merged,reduction = "umap",label = TRUE,pt.size = 1.5)
DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)

DimPlot(merged,reduction = "umap",label = TRUE,group.by="donor",pt.size = 0.7)
DimPlot(merged,reduction = "umap",label = TRUE,group.by="orig.ident",pt.size = 0.7)
DimPlot(merged,reduction = "umap",label = TRUE,split.by="donor",pt.size = 0.7)
DimPlot(merged,reduction = "umap",label = TRUE,split.by="orig.ident",pt.size = 0.7)

DefaultAssay(merged) <- "RNA"
dim(merged)
markers <- FindAllMarkers(
  merged,
  only.pos = TRUE, # 如果是TRUE,则只显示上调的差异基因
  min.pct = 0.25, # 如果是做GSEA,这里改成0
  logfc.threshold = 0.25 # 如果是做GSEA,这里改成0
)
head(markers)
write.table(markers,"all_cluster_markers0.3_unanno.csv")
saveRDS(merged,file = "seurat_dim30_0.3_unanno.rds")

new.cluster.ids <- c("Smooth muscle cell", "Macrophage","Fibroblast", 
                     "Fibroblast", "Smooth muscle cell", "T cell",
                     "Endothelial","Smooth muscle cell", 
                     "NK cell","Macrophage", "Myofibroblast",
                     "T cell", "Endothelial", "Mast cell", 
                     "Fibroblast","Macrophage","Endothelial", "Oligodendrocyte"
)
names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
merged@meta.data$cell_type <-Idents(merged)
table(merged$cell_type)

DimPlot(merged, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 0.5)
saveRDS(merged,file = "seurat_dim30_0.3_anno.rds")

##
control <- readRDS("./seurat_dim30_0.3_anno.rds")
control$sample <- control$orig.ident
control$group <- "control"

case <- readRDS("./seurat_integration_anno2.rds")
case$sample <- case$orig.ident
case$sample <- case$donor
case$donor[which(case$donor == c("rac1-1")) ] <- "rac1"
case$donor[which(case$donor == c("rac1-2")) ] <- "rac1"
case$donor[which(case$donor == c("rac2-1")) ] <- "rac2"
case$donor[which(case$donor == c("rac2-2")) ] <- "rac2"
case$donor[which(case$donor == c("rac3-1")) ] <- "rac3"
case$donor[which(case$donor == c("rac3-2")) ] <- "rac3"
case$donor[which(case$donor == c("rac3-3")) ] <- "rac3"
case$donor[which(case$donor == c("rac4-1")) ] <- "rac4"
case$group <- "case"

merge <- merge(control,case)
merge
saveRDS(merge,"mergeall.rds")
head(merge)


##
seurat_obj <- readRDS("./mergeall.rds")
#### 质量控制
seurat_obj <- PercentageFeatureSet(seurat_obj,
                                   pattern = "^MT-",
                                   col.name = "percent.mt")

seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
head(seurat_obj$percent.rb)
Idents(seurat_obj) <- "donor"

VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rb"),
  pt.size = 0,ncol = 4
)
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rb"),
  pt.size = 0.5,ncol = 4
)

seurat_obj <- subset(seurat_obj, subset =
                       nFeature_RNA >= 200 &
                       nFeature_RNA <= 3500 &
                       percent.mt <= 10)
dim(seurat_obj)

VlnPlot(seurat_obj,features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0,)
VlnPlot(seurat_obj,features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.rb"),ncol = 4, pt.size = 0.5)

seurat_obj <- FindVariableFeatures(seurat_obj,
                       selection.method = "vst",
                       nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj,vars.to.regress = c("nFeature_RNA","percent.mt","percent.rb"))  
seurat_obj <- RunPCA(seurat_obj)

DimPlot(seurat_obj,reduction = "pca", group.by = "donor")
sce <- seurat_obj

library(harmony)
system.time({sce <- RunHarmony(sce, group.by.vars = c("sample"))})
harmony_embeddings <- Embeddings(sce, 'harmony')
harmony_embeddings[1:5, 1:5]
DimPlot(object = sce, reduction = "harmony", pt.size = .1, group.by = "donor")
VlnPlot(object = sce, features = "harmony_1", group.by = "donor", pt.size = .1)

ElbowPlot(seurat_obj, ndims = 50)
ndim = 30 
sce <- RunUMAP(sce,  dims = 1:ndim,reduction = "harmony")
sce <- RunTSNE(sce,  dims = 1:ndim,reduction = "harmony")
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:ndim)
sce <- FindClusters(sce, resolution = 0.6)

DimPlot(sce,reduction = "umap",label = T,pt.size = 0.2)
DimPlot(sce,reduction = "tsne",label = T,pt.size = 0.2)
DimPlot(sce,group.by = "donor",reduction = "umap",label = T,pt.size = 0.2)
DimPlot(sce,split.by ="donor",reduction = "umap",label = T,pt.size = 0.2)
saveRDS(sce,"./1-unanno.rds")

##
seurat <- readRDS("1-unanno.rds")
dim(seurat)
new.cluster.ids <- c(
  "Smooth muscle cell", # 00  
  "Macrophage",# 01  
  "Fibroblast", # 02  
  "Fibroblast", # 03  
  "T cell", # 04 
  "Smooth muscle cell",# 05 
  "Endothelial",# 06
  "Fibroblast", #07 
  "Smooth muscle cell",#08 
  "Macrophage",#09
  "NK cell",#10
  "Fibroblast",#11 
  "Macrophage",#12
  "Endothelial",#13 
  "Fibroblast",#14
  "Myofibroblast",#15 
  "Mast cell",#16 
  "B cell",#17
  "Plasma cell",#18 
  "Macrophage",#19
  "Neuron",#20
  "Macrophage",#21
  "Smooth muscle cell",#22
  "Macrophage",#23
  "Fibroblast",#24
  "Pericyte",#25
  "Endothelial"#26
)
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj,
                           new.cluster.ids)

seurat.obj$cell_type <- Idents(seurat.obj)

saveRDS(seurat.obj, file = "2_anno.rds")


