mamba create -n auto-cell-type-annotation -c conda-forge -c bioconda -c defaults r-base r-devtools r-biocmanager -y

mamba install bioconductor-celldex bioconductor-singlecellexperiment r-seurat bioconductor-singler -y

mamba activate auto-cell-type-annotation