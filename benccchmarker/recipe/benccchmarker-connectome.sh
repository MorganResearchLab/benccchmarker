# Requires Seurat v3 and r-base=3.6.3

mamba create -n benccchmarker-connectome -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.0 r-devtools r-biocmanager r-ggplot2 r-seurat=3.2.3 cmake r-circlize bioconductor-complexheatmap r-igraph r-leiden pandas r-cpp11 -y

mamba activate benccchmarker-connectome

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"


Rscript -e "install.packages('igraph')"
Rscript -e "install.packages('plotrix')"
Rscript -e "devtools::install_github('msraredon/Connectome', ref = 'master')"
Rscript -e "install.packages('remotes'); remotes::install_version('SeuratObject', version = '4.0.0'); remotes::install_version('Seurat', version = '4.0.0'); remotes::install_version('uwot', version = '0.1.10'); remotes::install_version('sctransform', version = '0.3.2')"
