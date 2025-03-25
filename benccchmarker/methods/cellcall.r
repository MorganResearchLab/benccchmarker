library(cellcall)
library(Seurat)
library(yaml)
library(tibble)
library(tidyr)
library(dplyr)
library(peakRAM)

start_time <- Sys.time()

mem <- peakRAM({
    # Load argument from command line
    args <- commandArgs(trailingOnly = TRUE)

    # Load the parameter config file
    config <- yaml::read_yaml(args[3])
    
    seurat_object <- readRDS(args[1])

    counts <- as.matrix(seurat_object[["RNA"]]$data)
    colnames(counts) <- gsub("_", "", colnames(counts))
    colnames(counts) <- paste(colnames(counts), gsub("-", " ", seurat_object$cell_type), sep = "_")

    # Remove - from the cell id change it to "" # nolint
    colnames(counts) <- gsub("-", "", colnames(counts))

    counts <- as.data.frame(counts)

    end_preprocessing_time <- Sys.time()

    if (args[2] == "hsapiens") {
        org <- "Homo sapiens"
    } else {
        org <- "Mus musculus"
    }

    mt <- CreateNichConObject(
        data = counts,
        min.feature = config$params$CreateNichConObject$min.feature,
        names.field = config$params$CreateNichConObject$names.field,
        names.delim = config$params$CreateNichConObject$names.delim,
        source = config$params$CreateNichConObject$source,
        scale.factor = config$params$CreateNichConObject$scale.factor,
        Org = org,
        project = config$params$CreateNichConObject$project
    ) # nolint: indentation_linter.

    mt <- TransCommuProfile(
        object = mt,
        pValueCor = config$params$TransCommuProfile$pValueCor,
        CorValue = config$params$TransCommuProfile$CorValue,
        topTargetCor = config$params$TransCommuProfile$topTargetCor,
        p.adjust = config$params$TransCommuProfile$p.adjust,
        use.type = config$params$TransCommuProfile$use.type,
        probs = config$params$TransCommuProfile$probs,
        method = config$params$TransCommuProfile$method,
        Org = org,
        IS_core = TRUE
    )

    end_inference_time <- Sys.time()

    result <- mt@data$expr_l_r_log2_scale

    result <- as.data.frame(result)
    result <- rownames_to_column(result, var = "LR")
    result <- pivot_longer(result, cols = -LR, names_to = "sr", values_to = "value")
    result <- result[which(result$value > 0), ]
    result <- separate(data = result, col = sr, into = c("source", "target"), sep = "-")
    result <- separate(data = result, col = LR, into = c("ligand", "receptor"), sep = "-")
    colnames(result)[5] <- "LRscore"
    result <- result[which(result$source != result$target), ]
    result$ligand <- gsub(",", "&", result$ligand)
    result$receptor <- gsub(",", "&", result$receptor)
    result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = "---")

    result$label <- tolower(result$label)
    result <- distinct(result, label, .keep_all = TRUE)

    write.csv(result, file = args[4], row.names = FALSE)
})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time=preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)
