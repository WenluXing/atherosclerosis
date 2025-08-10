setwd("/home/xingwl/")
set.seed(717)
#### Load packages ----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

#### Configuration ----
# Seurat object rds path?
seurat_obj_path <- "/home/xingwl/practice5-Atherosis/noncoding/anno/harmony_anno_20_0.8.rds"
# Ready to analyze cell type?
selected_cell_type <- "Endothelial"   #758
# Number of nearest neighbors to aggregate?
# k <- 10  #668
k <- 20  #225
# k <- 15  #551
# Output path?
outdir <- "./practice5-Atherosis/noncoding/results-V3-metacell"

#### Preprocessing ----
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "data")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj <- seurat_obj_path %>% readRDS()
sub_seurat_obj <- seurat_obj %>%
  subset(cell_type == selected_cell_type)

#### Construct metacells ----
sub_seurat_obj <- sub_seurat_obj %>%
  cat_construct_metacells(k = k, name = selected_cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(sub_seurat_obj)
table(metacell_obj$donor)
#### Save ----
sub_seurat_obj %>% saveRDS(file.path(outdir, "seurat_obj.rds"))
head(sub_seurat_obj)


