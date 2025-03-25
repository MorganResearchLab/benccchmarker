library(Seurat)
library(Matrix)
library(parallel)
library(scMLnet)
library(tidyr)
library(peakRAM)
library(Seurat)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({

    # Load the parameter config file
    config <- yaml::read_yaml(args[3])
    seurat_object <- readRDS(args[1])

    gene_expression_matrix <- seurat_object[["RNA"]]$counts
    GCMat <- as(gene_expression_matrix, "dgCMatrix")

    BarCluTable <- data.frame(Barcode = colnames(seurat_object), Cluster = seurat_object$cell_type)
    cell_types <- sort(unique(BarCluTable$Cluster))

    sender_ct <- cell_types
    receiver_ct <- cell_types

    result <- list()

    for (LigClu in sender_ct) {
        for (RecClu in receiver_ct[which(receiver_ct != LigClu)]) {
            tryCatch({
                netList <- RunMLnet(
                    GCMat, 
                    BarCluFile = "BarCluTable.txt", 
                    RecClu = RecClu, 
                    LigClu = LigClu,
                    logfc = 0.15,
                    LigRecLib = "https://raw.githubusercontent.com/SunXQlab/scMLnet/master/database/LigRec.txt", 
                    TFTarLib = "https://raw.githubusercontent.com/SunXQlab/scMLnet/master/database/TFTargetGene.txt",
                    RecTFLib = "https://raw.githubusercontent.com/SunXQlab/scMLnet/master/database/RecTF.txt"
                )
                list.names <- paste(LigClu, RecClu, sep = "_")
                result[[list.names]] <- netList
            }, error = function(e) {
                message(sprintf("Error in RunMLnet for LigClu: %s, RecClu: %s", LigClu, RecClu))
                message(e)
            })
        }
    }

    result[which(is.na(result))] <- NULL

    result <- lapply(result, function(res){
        res <- res$LigRec
        res
    })
    result <- do.call(rbind, result)

    # Convert result to data frame
    result <- as.data.frame(result)
})