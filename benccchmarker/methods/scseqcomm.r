library(scSeqComm)
library(Seurat)
library(peakRAM)

args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({
    config <- yaml::read_yaml(args[3])

    seurat_object <- readRDS(args[1])

    gene_expression_matrix <- as.matrix(seurat_object[["RNA"]]$counts)
    cell_metadata <- as.data.frame(seurat_object$cell_type)
    cell_metadata <- data.frame(Cell_ID = rownames(cell_metadata), Cluster_ID = cell_metadata[,1])
    colnames(cell_metadata)[2] <- "Cluster_ID"

    LR_db <- LR_pairs_Kumar_2018
    TF_TG_db <- TF_TG_TRRUSTv2_HTRIdb_RegNetwork_High
    TF_PPR <- TF_PPR_KEGG_human
    
    end_preprocessing_time <- Sys.time()

    print(config$params$scSeqComm_analyze)

    comm_res <- scSeqComm_analyze(
        gene_expr = gene_expression_matrix,
        cell_metadata = cell_metadata,
        LR_pairs_DB = LR_db,
        TF_reg_DB = TF_TG_db,
        R_TF_association = TF_PPR,
        N_cores = config$params$N_cores,
        # inter_signaling = config$params$inter_signaling,
        # intra_signaling = config$params$intra_signaling,
        # DEmethod = config$params$DEmethod,
        # only.pos= config$params$only.pos,
        # aggregation_method = config$params$aggregation_method,
        # inter_scores = config$params$inter_scores,
        # sampling = config$params$sampling,
        # Nrep = config$params$Nrep,
        # min_cells = config$params$min_cells,
        # count_thr = config$params$count_thr,
    )

    inter_intra_scores <- comm_res$comm_results

    inter_max_intra_scores <- scSeqComm_summaryze_S_intra(inter_intra_scores)

    result <- scSeqComm_select(inter_max_intra_scores, 
                                    S_inter = config$params$result_filter$S_inter, 
                                    S_intra = config$params$result_filter$S_intra,)

    end_inference_time <- Sys.time()

    colnames(result) <- tolower(colnames(result))

    # Rename column cluster_L in result to "source" and cluster_R to "target"
    colnames(result)[colnames(result) == "cluster_l"] <- "source"
    colnames(result)[colnames(result) == "cluster_r"] <- "target"

    result$label <- paste(result$source, result$ligand, result$target, result$receptor, sep = '---')

    # Save the memory and time used
    write.csv(result, args[4], row.names = FALSE)
  
})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time = preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)