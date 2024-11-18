library(Seurat)
library(SeuratDisk)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

seurat_obj <- readRDS(args[1])

seurat_object[["RNA3"]] <- as(object = seurat_object[["RNA"]], Class = "Assay")
DefaultAssay(seurat_object) <- "RNA3"
seurat_object[["RNA"]] <- NULL
seurat_object <- RenameAssays(seurat_object, RNA3 = 'RNA')

save_file <- args[2]
Convert(save_file, dest = "h5ad")