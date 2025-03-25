mamba create -n benccchmarker-ccinx -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.4.2 r-devtools bioconductor-singlecellexperiment bioconductor-summarizedexperiment r-seurat r-dplyr -y

mamba activate benccchmarker-ccinx

Rscript -e "devtools::install_github('BaderLab/scClustViz')"
Rscript -e "devtools::install_github('immunogenomics/presto')"
Rscript -e "devtools::install_github('BaderLab/CCInx')"