library(SingleCellSignalR)
library(Seurat)
library(peakRAM)

start_time <- Sys.time()

mem <- peakRAM({
    args <- commandArgs(trailingOnly = TRUE)

    config <- yaml::read_yaml(args[3])

    seurat_object <- readRDS(args[1])

    i <- 1
    for (ct in unique(seurat_object$cell_type)) {
        seurat_object[[]][which(seurat_object$cell_type == ct), "cluster"] <- i
        i <- i+1
    }

    data <- seurat_object[["RNA"]]$counts

    if (args[2] == "hsapiens") {
        species <- "homo sapiens"
    } else {
        species <- "mus musculus"
    }

    cluster_names <- as.character(unique(seurat_object$cell_type))

    signal <- cell_signaling(
        data = data,
        genes = rownames(data),
        int.type = config$params$cell_signaling$int.type,
        species = species,
        cluster = seurat_object$cluster,
        c.names = cluster_names,
        s.score = config$params$cell_signaling$s.score,
        logFC = log2(config$params$cell_signaling$logFC),
        write = config$params$cell_signaling$write,
        verbose = config$params$cell_signaling$verbose,
    )

    intercellular_network <- tryCatch({
        inter_network(
            data = data,
            signal = signal,
            genes = rownames(data),
            cluster = seurat_object$cluster,
            c.names = cluster_names,
            species = species,
            write = FALSE
        )
    }, error = function(e) {
        message("Error when generating paracrine intercellular network", e$message)
        NULL # Return NULL if there's an error
    })

    result <- intercellular_network$`full-network`
    result$ligand <- gsub("\\.", "---", result$ligand)
    result$receptor <- gsub("\\.", "---", result$receptor)

    result$label <- paste(result$ligand, result$receptor, sep = "---")

    result$label <- tolower(result$label)

    # Split "ligand" columns into "source" and "target" columns by "---"
    result$source <- sapply(strsplit(result$ligand, "---"), function(x) x[1])
    result$ligand <- sapply(strsplit(result$ligand, "---"), function(x) x[2])
    result$target <- sapply(strsplit(result$receptor, "---"), function(x) x[1])
    result$receptor <- sapply(strsplit(result$receptor, "---"), function(x) x[2])

    write.csv(result, file=args[4])
})

end_time <- Sys.time()
elapsed_time <- end_time - start_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(Peak_RAM_Used_MiB = mem$Peak_RAM_Used_MiB, Elapsed_Time = elapsed_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)

