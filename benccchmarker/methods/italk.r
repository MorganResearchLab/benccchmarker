library(iTALK)
library(Seurat)
library(peakRAM)


# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({
  # Load the parameter config file
  config <- yaml::read_yaml(args[3])

  seurat_object <- readRDS(args[1])

  matrix <- as.data.frame(t(as.matrix(seurat_object[["RNA"]]$counts)))
  matrix$cell_type <- seurat_object$cell_type

  highly_exprs_genes <- rawParse(
    matrix,
    top_genes = config$params$rawParse$top_genes, 
    stats = config$params$rawParse$stats
  )

	end_preprocessing_time <- Sys.time()

  comm_list <- config$communication_categories

  result <- NULL
  for (comm.type in comm_list) {
    res.tmp <- FindLR(
      highly_exprs_genes, 
      datatype = config$params$FindLR$datatype, 
      comm_type = comm.type
    )
    res.tmp <- res.tmp[order(res.tmp$cell_from_mean_exprs * res.tmp$cell_to_mean_exprs, decreasing = T), ]
    result <- rbind(result, res.tmp)
  }

	end_inference_time <- Sys.time()

  result$source <- result$cell_from
  result$target <- result$cell_to

  result$cell_from <- NULL
  result$cell_to <- NULL

  result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = "---")

  result <- result[!duplicated(result$label), ]

  write.csv(result, file = args[4], row.names = FALSE)
})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time = preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)
