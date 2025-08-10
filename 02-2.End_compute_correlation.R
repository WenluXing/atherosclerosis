setwd("/home/xingwl/")
#### Load packages ----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

#### Configuration -----
# Seurat object rds path?
seurat_obj <- readRDS("/home/xingwl/practice5-Atherosis/noncoding/results-V3-metacell/Endothelial/k_20/data/seurat_obj.rds")

# Ready to analyze genes?
selected_genes <- c("PSMB1","ZFAS1","SNHG32","SNHG5","SNHG6","SNHG8","SNHG29","NEAT1","MALAT1","GAS5",
                    "SNHG7","OTUD6B-AS1","PHF1","NORAD","RGS5","FTX","SNHG9","SNHG16","CYTOR",
                    "EPB41L4A-AS1",
                    
                    "LINC00963","PCAT19","MEG3","LINC02802","MIR99AHG",
                    "MAGI2-AS3","MIR22HG","KCNQ1OT1","CARMN","LINC01615",
                    "CRNDE","SERTAD4-AS1","MIR100HG","DANCR","CD27-AS1","PCED1B-AS1",
                    "LINC01871","PRKCQ-AS1","LINC00623","LINC-PINT","LINC00924",
                    "LINC01781","LINC00926","ILF3-DT","FAM215B","LINC02362","CHL1-AS2",
                    "FOXD3-AS1","NSMCE1-DT",
                    
                    "FTO","METTL3","METTL14","RBM15","RBM15B","WTAP","KIAA1429","CBLL1",
                    "ZC3H13","ALKBH5","YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1","IGF2BP2",
                    "IGF2BP3","HNRNPA2B1","HNRNPC","FMR1","LRPPRC","ELAVL1"
) 

# Ready to analyze cell type?
selected_cell_type <- "Endothelial"
#seurat.obj <- subset(seurat.obj,subset = cell_type %in% selected_cell_type)
# Number of nearest neighbors to aggregate?
k <- 20 #225
# k <- 15  #551
# Output path?
outdir <- "./practice5-Atherosis/noncoding/results-V4-correlation"

#### Preprocessing ----
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "correlation")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}

#seurat_obj <- seurat_obj_path %>% readRDS()
metacell_obj <- seurat_obj %>% hdWGCNA::GetMetacellObject()
metacell_obj %>% dim()
cell_types <- names(table(metacell_obj$cell_type))
selected_genes <-
  selected_genes[selected_genes %in% rownames(metacell_obj)]

#### Correlations between genes ----
all_genes <- rownames(metacell_obj)
genes_cor_res <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = all_genes)

# Save
write_csv(genes_cor_res,
          file = file.path(outdir,
                           "gene.cor.csv"))

#### Correlations between genes and pathway scores ----
category <- "H"   #关注啥通路，写啥通路
pathway_cor_res <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = category,
                  outdir = outdir,  #加上这一行会出一个打分的（后面画散点图要用到横坐标基因，纵坐标通路）
                  # species = "mouse"
  )   #默认物种是人 #默认打分是aucell比较快，score_method = "aucell"/"gsva"/"ssgsea"/"addmodulescore"

# Save
write_csv(pathway_cor_res,
          file = file.path(outdir,
                           str_c(category, ".cor.csv")))

# #### Correlations between genes and pathway scores ----
category <- "GO:BP"
pathway_cor_res1 <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = category,
                  outdir = outdir,
                  # species = "mouse"
  )

# Save
write_csv(pathway_cor_res1,
          file = file.path(outdir,
                           str_c(category, ".cor.csv")))
# 
#### Correlations between genes and pathway scores ----
category <- "CP:KEGG"
pathway_cor_res2 <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = category,
                  outdir = outdir,
                  # species = "mouse"
  )

# Save
write_csv(pathway_cor_res2,
          file = file.path(outdir,
                           str_c(category, ".cor.csv")))

