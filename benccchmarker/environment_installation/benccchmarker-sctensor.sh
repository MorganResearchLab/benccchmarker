# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

mamba create -n benccchmarker-sctensor -c conda-forge -c bioconda -c defaults -c python=3.10 r-base=4.3.3 r-devtools r-biocmanager r-ggplot2 r-seurat bioconductor cmake bioconductor-schex -y

mamba activate benccchmarker-sctensor

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

Rscript -e 'install.packages(c("RSQLite", "igraph", "plotly", "nnTensor", "rTensor", "abind", "plotrix", "heatmaply", "tagcloud", "rmarkdown", "knitr", "outliers", "crayon", "checkmate", "testthat", "Seurat", "BiocManager", "concaveman"), repos="http://cran.r-project.org")'

Rscript -e 'BiocManager::install(c("S4Vectors", "reactome.db", "AnnotationDbi", "SummarizedExperiment", "SingleCellExperiment", "BiocStyle", "biomaRt", "MeSHDbi", "Category", "meshr", "GOstats", "ReactomePA", "DOSE", "LRBase.Hsa.eg.db", "MeSH.Hsa.eg.db", "LRBase.Mmu.eg.db", "MeSH.Mmu.eg.db", "LRBaseDbi", "Homo.sapiens", "schex"), suppressUpdates=TRUE)'

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src

Rscript -e 'devtools::install_github("rikenbit/scTensor")'