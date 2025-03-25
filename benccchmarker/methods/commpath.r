library(CommPath)
library(peakRAM)
library(Seurat)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({

    # Load the parameter config file
    config <- yaml::read_yaml(args[3])
    seurat_object <- readRDS(args[1])

    species <- args[2]

    # Raise an error if the species is not in ('hsapies', 'mmusculus', 'rnorvegicus')
    if (!(species %in% c('hsapiens', 'mmusculus', 'rnorvegicus'))) {
        stop("The species must be one of 'hsapiens', 'mmusculus', 'rnorvegicus'")
    }

    commpath_object <- createCommPath(
        expr.mat = seurat_object[["RNA"]]$data, 
		cell.info = seurat_object$cell_type, 
		species = species
    )

    end_preprocessing_time <- Sys.time()

    commpath_object <- findLRmarker(
        object = commpath_object, 
        method = config$params$findLRmarker$method,
    )
    commpath_object <- findLRpairs(
        object = commpath_object,
		logFC.thre = config$params$findLRpairs$logFC.thre, 
		p.thre = config$params$findLRpairs$p.thre
    )

    end_inference_time <- Sys.time()
    
    result <- commpath_object@interact[['InteractGene']]

    result$source <- result$cell.from
    result$target <- result$cell.to
    result$cell.from <- NULL
    result$cell.to <- NULL

    result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = '---')

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