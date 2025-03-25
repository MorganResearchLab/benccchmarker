mamba create -n benccchmarker-commpath -c conda-forge -c bioconda -c defaults python=3.9 r-base=4.3.1 pip r-matrix r-circlize r-ggplot2 r-dplyr r-reshape2 bioconductor-gsva r-devtools r-seurat r-seuratobject -y 

mamba activate benccchmarker-commpath

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"


Rscript -e 'devtools::install_github("yingyonghui/CommPath")'