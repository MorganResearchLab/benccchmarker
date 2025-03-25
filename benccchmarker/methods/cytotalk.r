library(CytoTalk)
library(Seurat)
library(yaml)
library(peakRAM)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({
    args <- commandArgs(trailingOnly = TRUE)

    config <- yaml::read_yaml(args[3])
    seurat_object <- readRDS(args[1])

    counts <- as.matrix(seurat_object[["RNA"]]$data)

    metadata <- data.frame(colnames(seurat_object), seurat_object$cell_type)
    colnames(metadata) <- c("barcode_sample", "cell_type")

    main_dir <- dirname(args[4])

    if (!dir.exists(file.path(main_dir, "tmp/cytotalk"))) {
        dir.create(file.path(main_dir, "tmp/cytotalk"))
    }

    tmp_count_path <- file.path(main_dir, "tmp/cytotalk/counts.csv")
    tmp_metadata_path <- file.path(main_dir, "tmp/cytotalk/metadata.csv")

    write.csv(counts, file = tmp_count_path, quote = FALSE)
    write.csv(metadata, file = tmp_metadata_path, quote = FALSE, row.names = FALSE)

    cytotalk_data_object <- read_matrix_with_meta(tmp_count_path, tmp_metadata_path)

    # Get unique cell types from cytotalk_data_object$cell_types
    cell_types <- unique(cytotalk_data_object$cell_types)

    end_preprocessing_time <- Sys.time()

    all_inference <- data.frame()

    for (type_a in cell_types) {
        for (type_b in cell_types) {
            tryCatch({
                results <- run_cytotalk(
                    cytotalk_data_object, 
                    type_a, 
                    type_b,
                    cutoff_a = config$params$run_cytotalk$cutoff_a,
                    cutoff_b = config$params$run_cytotalk$cutoff_b,
                    beta_max = config$params$run_cytotalk$beta_max,
                    omega_min = config$params$run_cytotalk$omega_min,
                    omega_max = config$params$run_cytotalk$omega_max,
                    depth = config$params$run_cytotalk$depth,
                    ntrial = config$params$run_cytotalk$ntrial,
                    cores = config$params$run_cytotalk$cores
                )

                inference <- results$pcst$final_network

                inference$label <- paste(inference$node1_type, inference$node1, inference$node2_type, inference$node2, sep = '---')

                all_inference <- rbind(all_inference, inference)
            }, error = function(e) {
                message("Error when inferencing ", type_a, " and ", type_b)
            })

        }
    }

    end_inference_time <- Sys.time()

    # Rename all_inference column node_1 to ligand, node_2 to receptor, node_1_type to source, node_2_type to target
    colnames(all_inference) <- c("ligand", "receptor", "source", "target", "node1_prize","node2_prize","node1_pem","node2_pem","is_ct_edge","cost","label")

    write.csv(all_inference, file = args[4], row.names = FALSE)

})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time=preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)



