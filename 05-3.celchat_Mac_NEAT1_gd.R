library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)

source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
seurat_obj <-
  readRDS("/home/xingwl/practice5-Atherosis/noncoding/anno/harmony_anno_20_0.8.rds")
#提取其余细胞
seurat_obj0 <- subset(seurat_obj,idents = c("Endothelial","Fibroblast","Fibromyocyte",
                                            "T cell","SMC","Pericyte","Plasma cell",
                                            "B cell","Neuron","NK cell","Mast cell"))#
#提取巨噬细胞高低表达分组
seurat_obj1 <- subset(seurat_obj,idents = c("Macrophage"))
selected_genes <- c("NEAT1")
expr <- seurat_obj1 %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj1[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj1[[str_c(selected_genes, "_group")]]))
table(Idents(seurat_obj1))
#将巨噬细胞的高低组提取出来
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
# cellchat输出文件路径
output.dir <-
  paste0("/home/xingwl/practice5-Atherosis/noncoding/results-V12-cellchat-Mac-NEAT1-gd/Mac_NEAT1_gd/") # , group, "/"must have "/"
dir.create(output.dir, recursive = T)

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

pdf(file = paste0(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()

pdf(file = paste0(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
  
)
dev.off()

pdf(file = paste0(output.dir, "netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
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
saveRDS(cellchat, file = paste0(output.dir, "cellchat.rds"))
cellchat <- readRDS("/home/xingwl/practice5-Atherosis/noncoding/results-V12-cellchat-Mac-NEAT1-gd/Mac_NEAT1_gd/cellchat.rds")



