library(CCInx)
library(SingleCellExperiment)
library(Seurat)
library(dplyr)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({
    # Load the parameter config file
    config <- yaml::read_yaml(args[3])
    species <- args[2]
    seurat_object <- readRDS(args[1])

    sce_object <- as.SingleCellExperiment(seurat_object)

    colnames(colData(sce_object))[colnames(colData(sce_object)) == "cell_type"] <- "cellTypes"

    end_preprocessing_time <- Sys.time()

    gsl <- BuildGeneStatList(
        inD = sce_object,
        cl = colData(sce_object)$cellTypes,
        assayType = config$params$BuildGeneStatList$assayType,
        assaySlot = config$params$BuildGeneStatList$assaySlot,
        exponent = config$params$BuildGeneStatList$exponent,
        pseudoconfig = config$params$BuildGeneStatList$pseudoconfig
    )

    inx <- BuildCCInx(
        GeneStatList=gsl,
        Species=species
    )

    end_inference_time <- Sys.time()

    edges_with_types <- inx$edges %>%
        left_join(inx$nodes %>% select(node, proteinType, MeanNormGeneExpr), by = c("nodeA" = "node")) %>%
        left_join(inx$nodes %>% select(node, proteinType, MeanNormGeneExpr), by = c("nodeB" = "node"), suffix = c(".A", ".B"))

    filtered_edges <- edges_with_types %>%
    filter((proteinType.A == "Receptor" & proteinType.B == "Ligand") | 
            (proteinType.A == "Ligand" & proteinType.B == "Receptor"))

    filtered_edges <- filtered_edges %>%
        arrange(desc(edgeWeight))

    result <- filtered_edges %>%
        mutate(Ligand = ifelse(proteinType.A == "Ligand", sapply(strsplit(nodeA, "_"), `[`, 1), 
                                ifelse(proteinType.B == "Ligand", sapply(strsplit(nodeB, "_"), `[`, 1), NA)),
                Receptor = ifelse(proteinType.A == "Receptor", sapply(strsplit(nodeA, "_"), `[`, 1), 
                                ifelse(proteinType.B == "Receptor", sapply(strsplit(nodeB, "_"), `[`, 1), NA)),
                Source = ifelse(proteinType.A == "Ligand", sapply(strsplit(nodeA, "_"), `[`, 2), 
                                ifelse(proteinType.B == "Ligand", sapply(strsplit(nodeB, "_"), `[`, 2), NA)),
                Target = ifelse(proteinType.A == "Receptor", sapply(strsplit(nodeA, "_"), `[`, 2), 
                                ifelse(proteinType.B == "Receptor", sapply(strsplit(nodeB, "_"), `[`, 2), NA))) %>%
        select(Ligand, Receptor, Source, Target, edgeWeight) %>%
        distinct()

    write.csv(result, file = args[4], row.names = FALSE)
})

preprocessing_time <- end_preprocessing_time - start_time
inference_time <- end_inference_time - end_preprocessing_time

met_path <- gsub(".csv", "_met.csv", args[4])

# Create a data frame with the peak RAM usage and elapsed time
mem <- data.frame(peak_ram_used_mib = mem$Peak_RAM_Used_MiB, preprocessing_time = preprocessing_time, inference_time = inference_time)

# Write the peak RAM usage and elapsed time to a file
write.csv(mem, file = met_path, row.names = FALSE)