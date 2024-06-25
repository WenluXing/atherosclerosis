#### Load packages ----
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
source("~/share/custom_function.R")
source("~/share/custom_plot_function.R")
output.dir <- paste0("~/results-V24-Figure/Fig7/Mac_NEAT1_meta") 

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

seurat_obj0 <- merge(Bcell1,c(End1,Fib1,Fibm1,NK1,Pericyte1,Plasma1,SMC1,Tcell1,Neu1))
seurat_obj1 <- Mac1

#
selected_genes <- c("NEAT1")
expr <- seurat_obj1 %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
# [1] 2.574677
seurat_obj1[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj1[[str_c(selected_genes, "_group")]]))
table(Idents(seurat_obj1))

#
group1 <- "low"
group2 <- "high"
seurat_obj2 <- seurat_obj1[, seurat_obj1[["NEAT1_group"]] == group1]
table(seurat_obj2$cell_type)
seurat_obj2$cell_type <- "Mac_NEAT1_low"
seurat_obj3 <- seurat_obj1[, seurat_obj1[["NEAT1_group"]] == group2]
table(seurat_obj3$cell_type)
seurat_obj3$cell_type <- "Mac_NEAT1_high"

merged <- merge(seurat_obj0, c(seurat_obj2, seurat_obj3))
table(merged$cell_type)
Idents(merged) <- "cell_type"
table(Idents(merged))

# cellchat
expr <- merged@assays$RNA@data
table(Idents(merged))

data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(merged))
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

pdf(file.path(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()

pdf(file.path(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
  
)
dev.off()

pdf(file.path(output.dir, "netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()
saveRDS(cellchat, file.path(output.dir, "cellchat.rds"))

##
cellchat <- readRDS("~/results-V23-Figure/Fig6/Mac_NEAT1_meta/cellchat.rds")
pathways.show <- c("CXCL")

#
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL"))
p <- gg1 + gg2
ggsave(file.path(output.dir, "san_CXCL_celltype.pdf"),
       plot = p,
       height = 4,
       width = 9.5)

#
par(mfrow=c(1,1))
pdf(file.path(output.dir,"heatmap_CXCL_celltype.pdf"))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds",width = 4, height = 3, font.size = 10)
dev.off()

#
pdf(file.path(output.dir,"chord_CXCL_gene_pre.pdf"))
netVisual_chord_gene(cellchat, sources.use = c(11:12), targets.use = c(2), signaling = c("CXCL"),legend.pos.x = 6)
dev.off()
pdf(file.path(output.dir,"chord_CXCL_gene.pdf"))
netVisual_chord_gene(cellchat, sources.use = c(11:12), targets.use = c(2), signaling = c("CXCL"),show.legend = FALSE)
dev.off()

#
pairLR.pathways.show <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.pathways.show[2,] 
vertex.receiver = seq(1,4)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
pdf("circlr_CXCL2-ACKR1_celltype.pdf")

LR.show <- pairLR.pathways.show[3,] 
vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
pdf("circlr_CXCL3-ACKR1_celltype.pdf")

LR.show <- pairLR.pathways.show[4,] 
vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
pdf("circlr_CXCL8-ACKR1_celltype.pdf")

#
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL"))
pdf(file.path(output.dir,"dotplot_Mac-Endo-CXCL.pdf"), height = 6, width = 2)
netVisual_bubble(cellchat, sources.use = c(11:12), targets.use = c(2), pairLR.use = pairLR.use, remove.isolate = TRUE, angle.x = 45)
dev.off()

#
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MIF","SPP1","VEGF","ANGPT","ANGPTL",
                                                        "EGF","FGF","PDGF","ncWNT","LAMININ","VISFATIN","ITGB2"))
pdf(file.path(output.dir,"dotplot_Mac-Endo-way.pdf"), height = 5, width = 3.5)
netVisual_bubble(cellchat, sources.use = c(11:12), targets.use = c(2), pairLR.use = pairLR.use, remove.isolate = TRUE, angle.x = 45)
dev.off()

#
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","NEGR","VEGF","PROS","NOTCH",
                                                        "IL16","IL10","ITGB2"))
pdf(file.path(output.dir,"dotplot_Endo-Mac-way.pdf"), height = 4.5, width = 3)
netVisual_bubble(cellchat, sources.use = c(2), targets.use = c(11:12), pairLR.use = pairLR.use, remove.isolate = TRUE, angle.x = 45)
dev.off()
