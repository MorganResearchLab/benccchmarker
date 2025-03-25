library(Seurat)
library(tidyr)
library(yaml)
library(CiteFuse)
library(SingleCellExperiment)
library(peakRAM)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({
    # Load the parameter config file
    config <- yaml::read_yaml(args[3])

    seurat_object <- readRDS(args[1])
    sce_object <- as.SingleCellExperiment(seurat_object, assay = "RNA")

    sce_object <- normaliseExprs(sce = sce_object,
        altExp_name = "none",
        exprs_value = "logcounts",
        transform = config$params$ligandReceptorTest$RNA_exprs_value
    )

    lr_list <- read.csv(config$params$ligandReceptorTest$ligandReceptor_list, header = TRUE, sep = ",")

    lr_list <- as.matrix(lr_list)
    lr_list <- matrix(unlist(strsplit(as.character(lr_list), " \\| ")), ncol = 2, byrow = TRUE)

    end_preprocessing_time <- Sys.time()

    sce_object <- ligandReceptorTest(
        sce = sce_object,
        ligandReceptor_list = lr_list,
        cluster = sce_object$cell_type,
        RNA_exprs_value = config$params$ligandReceptorTest$RNA_exprs_value,
        use_alt_exp = config$params$ligandReceptorTest$use_alt_exp,
        num_permute = config$params$ligandReceptorTest$num_permute,
    )

    end_inference_time <- Sys.time()

    cell_types <- levels(sce_object$cell_type)

    result <- as.data.frame(metadata(sce_object)$LRanalysis_results_RNA$LRanalysis_pvalue)
    result$lr <- rownames(result)
    write.csv(result, file = args[4], row.names = FALSE)
    
    result <- result %>%
        pivot_longer(
            cols = matches("\\|"),
            names_to = "st", 
            values_to = "p_value"
        )
    
    result <- result[result$p_value < config$p_value, ]

    lr_split <- do.call(rbind, strsplit(result$lr, "\\|"))
    result$ligand <- lr_split[, 1]
    result$receptor <- lr_split[, 2]

    st_split <- do.call(rbind, strsplit(result$st, "\\|"))
    result$source <- st_split[, 1]
    result$target <- st_split[, 2]

    result$lr <- NULL
    result$st <- NULL

    result$source <- cell_types[as.numeric(result$source)]
    result$target <- cell_types[as.numeric(result$target)]

    result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = '---')

    write.csv(result, file = args[4], row.names = FALSE)
})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time=preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)

