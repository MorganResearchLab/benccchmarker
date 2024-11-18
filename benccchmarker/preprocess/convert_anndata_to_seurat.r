library(Seurat)
library(anndata)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

ad <- read_h5ad(args[1])

counts <- t(as.matrix(ad$X))

seurat_obj <- CreateSeuratObject(
    counts,
    project = "simulated_data",
    meta.data = ad$obs
)

seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
)

# Save seurat_obj to rds
saveRDS(seurat_obj, args[2])
