library(Seurat)
library(celltalker)
library(reticulate)
library(tidyr)
library(peakRAM)


args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({

    config <- yaml::read_yaml(args[3])

    seurat_object <- readRDS(args[1])

    end_preprocessing_time <- Sys.time()

    result <- celltalk(
        seurat_object,
        metadata_grouping = "cell_type",
        ligand_receptor_pairs = ramilowski_pairs,
        number_cells_required = config$params$celltalk$number_cells_required,
        min_expression = config$params$celltalk$min_expression,
        max_expression = config$params$celltalk$max_expression,
        scramble_times = config$params$celltalk$scramble_times
    )
    result <- result[result$p_val < 0.05, ]

    end_inference_time <- Sys.time()

    # Split interaction_pairs column into source and target columns
    result <- separate(data = result, col = interaction_pairs, into = c("source", "target"), sep = "_")

    # Split interaction column into ligand and receptor columns
    result <- separate(data = result, col = interaction, into = c("ligand", "receptor"), sep = "_")

    result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = "---")

    result$label <- tolower(result$label)

    write.csv(result, file = args[4], row.names = FALSE)
})


preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time=preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)
