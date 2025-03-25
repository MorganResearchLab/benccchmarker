library(CellChat)
library(patchwork)
library(yaml)
library(Seurat)
options(stringsAsFactors = FALSE)
library(peakRAM)

# Load argument from command line
args <- commandArgs(trailingOnly = TRUE)
start_time <- Sys.time()

mem <- peakRAM({

    # Load the parameter config file
    config <- yaml::read_yaml(args[3])

    seurat_object <- readRDS(args[1])

    end_preprocessing_time <- Sys.time()

    cellchat <- createCellChat(
        object = seurat_object, 
        group.by = "cell_type"
    )

    if (args[2] == "hsapiens") {
        CellChatDB <- CellChatDB.human
    } else {
        CellChatDB <- CellChatDB.mouse
    }

    CellChatDB.use <- CellChatDB
    cellchat@DB <- CellChatDB.use

    cellchat <- subsetData(cellchat)

    future::plan(
        "multisession", 
        workers = config$general$workers
    )

    cellchat <- updateCellChat(cellchat)
    cellchat <- identifyOverExpressedGenes(
        cellchat,
        group.by = config$params$identifyOverExpressedGenes$group.by,
        idents.use = config$params$identifyOverExpressedGenes$idents.use,
        invert = config$params$identifyOverExpressedGenes$invert,
        group.dataset = config$params$identifyOverExpressedGenes$group.dataset,
        pos.dataset = config$params$identifyOverExpressedGenes$pos.dataset,
        features.name = config$params$identifyOverExpressedGenes$features.name,
        only.pos = config$params$identifyOverExpressedGenes$only.pos,
        features = config$params$identifyOverExpressedGenes$features,
        return.object = config$params$identifyOverExpressedGenes$return.object,
        thresh.pc = config$params$identifyOverExpressedGenes$thresh.pc,
        thresh.fc = config$params$identifyOverExpressedGenes$thresh.fc,
        thresh.p = config$params$identifyOverExpressedGenes$thresh.p,
    )

    cellchat <- identifyOverExpressedInteractions(cellchat)

    cellchat <- computeCommunProb(
        cellchat,
        type="triMean"
    )

    cellchat <- filterCommunication(
        cellchat, 
        min.cells = config$params$filterCommunication$min.cells
    )

    end_inference_time <- Sys.time()

    result <- subsetCommunication(cellchat)

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