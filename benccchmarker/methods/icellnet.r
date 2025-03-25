library(BiocGenerics)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(jetset)
library(ggplot2)
library(dplyr)
library(icellnet)
library(gridExtra)
library(Seurat)
library(peakRAM)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({
  # Load the parameter config file
  config <- yaml::read_yaml(args[3])

  seurat_object <- readRDS(args[1])
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- ScaleData(seurat_object)

  database_path <- "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"

  database <- read.csv(
    database_path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = ""
  )
  database <- as.data.frame(database)

  database_coupled <- name.lr.couple(database, type = "Family")

  Idents(seurat_object) <- seurat_object$cell_type

  preprocessed_input <- sc.data.cleaning(
    object = seurat_object,
    db = database,
    filter.perc = config$params$sc.data.cleaning$filter.perc,
    save_file = FALSE
  )

  rescaled_input <- as.data.frame(gene.scaling(
    preprocessed_input,
    n = config$params$gene.scaling$n_var,
    db = database
  ))

  pc_data <- rescaled_input

  pc_target <- data.frame(
    Class = colnames(pc_data)[-dim(rescaled_input)[2]],
    ID = colnames(pc_data)[-dim(rescaled_input)[2]],
    Cell_type = colnames(pc_data)[-dim(rescaled_input)[2]]
  )
  rownames(pc_target) <- colnames(pc_data)[-dim(rescaled_input)[2]]

  pc_ct <- colnames(pc_data)[-dim(rescaled_input)[2]]
  cc_ct <- colnames(pc_data)[-dim(rescaled_input)[2]]

  end_preprocessing_time <- Sys.time()

  results <- lapply(cc_ct, function(ct) {
    score.computation <- icellnet.score(
      direction = "out",
      PC.data = pc_data,
      CC.data = as.data.frame(rescaled_input[, ct],
        row.names = rownames(rescaled_input)
      ),
      PC.target = pc_target,
      PC = pc_ct[which(pc_ct != ct)],
      CC.type = "RNAseq",
      PC.type = "RNAseq",
      db = database
    )
    lr <- as.matrix(score.computation[[2]][apply(score.computation[[2]], 1, function(y) any(!is.na(y))), ])
    lr <- as.matrix(lr[which(rowSums(lr) > 0), ])
    lr
  })

  end_inference_time <- Sys.time()

  names(results) <- colnames(pc_data)[-dim(rescaled_input)[2]]

  result <- lapply(names(results), function(ct) {
    temp_result <- results[[ct]]
    temp_result <- as.data.frame(temp_result)
    temp_result <- tibble::rownames_to_column(temp_result, "ligand_receptor")

    temp_result <- temp_result %>% tidyr::pivot_longer(cols = -ligand_receptor, names_to = "target", values_to = "LRscore")

    temp_result$source <- ct
    temp_result <- tidyr::separate(data = temp_result, col = ligand_receptor, into = c("ligand", "receptor"), sep = " / ")

    temp_result <- temp_result %>%
      tidyr::separate_rows(ligand, sep = " \\+ ") %>%
      tidyr::separate_rows(receptor, sep = " \\+ ")
  })

  result <- do.call(rbind, result)

  result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = "---")

  write.csv(result, file = args[4], row.names = FALSE)
})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time = preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)








write.csv(result, file = "ICELLNET.csv", row.names = F)
