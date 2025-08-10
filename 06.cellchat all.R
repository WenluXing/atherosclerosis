#### Fig6A FigS6A ----
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
output.dir <- paste0("./results-V23-Figure/Fig6/all_meta") 

##
Bcell <- readRDS("./results-V3-metacell/B cell/k_20/data/seurat_obj.rds")
Bcell1 <- hdWGCNA::GetMetacellObject(Bcell)
End <- readRDS("./results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")
End1 <- hdWGCNA::GetMetacellObject(End)
Fib <- readRDS("./results-V3-metacell/Fibroblast/k_20/data/seurat_obj.rds")
Fib1 <- hdWGCNA::GetMetacellObject(Fib)
Fibm <- readRDS("./results-V3-metacell/Fibromyocyte/k_20/data/seurat_obj.rds")
Fibm1 <- hdWGCNA::GetMetacellObject(Fibm)
Mac <- readRDS("./results-V3-metacell/Macrophage/k_20/data/seurat_obj.rds")
Mac1 <- hdWGCNA::GetMetacellObject(Mac)
NK <- readRDS("./results-V3-metacell/NK cell/k_20/data/seurat_obj.rds")
NK1 <- hdWGCNA::GetMetacellObject(NK)
Pericyte <- readRDS("./results-V3-metacell/Pericyte/k_20/data/seurat_obj.rds")
Pericyte1 <- hdWGCNA::GetMetacellObject(Pericyte)
Plasma <- readRDS("./results-V3-metacell/Plasma cell/k_20/data/seurat_obj.rds")
Plasma1 <- hdWGCNA::GetMetacellObject(Plasma)
SMC <- readRDS("./results-V3-metacell/SMC/k_20/data/seurat_obj.rds")
SMC1 <- hdWGCNA::GetMetacellObject(SMC)
Tcell <- readRDS("./results-V3-metacell/T cell/k_20/data/seurat_obj.rds")
Tcell1 <- hdWGCNA::GetMetacellObject(Tcell)
Neu <- readRDS("./results-V3-metacell/Neuron/k_20/data/seurat_obj.rds")
Neu1 <- hdWGCNA::GetMetacellObject(Neu)

seurat_obj <- merge(Bcell1,c(End1,Fib1,Fibm1,NK1,Pericyte1,Plasma1,SMC1,Tcell1,Neu1,Mac1))
Idents(seurat_obj) <- seurat_obj$cell_type

expr <- seurat_obj@assays$RNA@data

data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(seurat_obj))
colnames(meta) <- "labels"

unique(meta$labels) # check the cell labels
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <-
  as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <-
  CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multisession", workers = 10) # do parallel
cellchat <-
  identifyOverExpressedGenes(cellchat) # take a short time
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) # take a short time

#### Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#### Extract the inferred cellular communication network as a data frame
#pass
#### Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#### Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

saveRDS(cellchat, file.path(output.dir, "cellchat.rds"))

#### Load data ----
cellchat <- readRDS("./results-V23-Figure/Fig6/all_meta/cellchat.rds")

#
par(mfrow = c(1,2), xpd=TRUE)
groupSize <- as.numeric(table(cellchat@idents))
pdf(file.path(output.dir, "FigS6A_circle_number.pdf"),height=6,width=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf(file.path(output.dir, "Fig6A_circle_weight.pdf"),height=6,width=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
