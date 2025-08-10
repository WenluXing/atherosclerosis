setwd("/home/xingwl/")
#### Load packages ----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

#### Configuration -----
# Seurat object rds path?
seurat_obj_path <-
  "/home/xingwl/practice5-Atherosis/noncoding/results-V3-metacell/SMC/k_20/data/seurat_obj.rds"
#阴阳分组使用原来的seurat对象，高低分组使用metacell之后的对象
#但是无论使用哪一种都是这样的路径
# Ready to analyze genes?
selected_genes <- c("NEAT1")  
# Ready to analyze cell type?
selected_cell_type <- "SMC"
# Number of nearest neighbors to aggregate?
k <- 20
# Output path?
outdir <- "/home/xingwl/practice5-Atherosis/noncoding/results-V4-correlation/"

# #### Preprocessing ----
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "NEAT1_deg")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

seurat_obj <- seurat_obj_path %>% readRDS()
seurat_obj
DimPlot(seurat_obj,group.by = "cell_type")

##------------------第2种 高低------------------------
options(future.globals.maxSize = 1000 * 1024^3)  #报错内存不够时添加此行
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
# ### Compute deg ----
# install.packages("~/R/future_1.30.0.tar.gz", repos = NULL, type = "source")
library(future)
deg_metacell_obj <- metacell_obj |> cat_deg(
  group_by = selected_genes,
  cell_types = selected_cell_type,
  min_pct = 0,
  mode = 'median',
  logfc_threshold = 0
)
table(metacell_obj@assays$RNA@counts['NEAT1', ] > median(metacell_obj@assays$RNA@counts['NEAT1', ]))
write_csv(deg_metacell_obj, file.path(outdir, "deg_metacell_obj.csv"))

#### Function enrich ----
enriched <- deg_metacell_obj |>
  filter(group == "NEAT1") |>
  cat_enrich()   
saveRDS(enriched, file.path(outdir, "NEAT1_enriched.rds"))

#### GSEA ----
gsea_res <- deg_metacell_obj |>
  filter(group == "NEAT1") |>
  cat_gsea(category = "H"
           # species = "mouse"
  )
saveRDS(gsea_res, file.path(outdir, "h_gsea.rds"))

gsea_res1 <- deg_metacell_obj |>
  filter(group == "NEAT1") |>
  cat_gsea(category = "C2"
           # species = "mouse"
  )
saveRDS(gsea_res1, file.path(outdir, "C2_gsea.rds"))

gsea_res2 <- deg_metacell_obj |>
  filter(group == "NEAT1") |>
  cat_gsea(category = "C5"
           # species = "mouse"
  )
saveRDS(gsea_res2, file.path(outdir, "C5_gsea.rds"))

