# Requires Seurat v3 and r-base=3.6.3

mamba create -n benccchmarker-connectome -c conda-forge -c bioconda -c defaults python=3.10 r-base=3.6.3 -y

mamba activate benccchmarker-connectome

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install -y r-devtools r-biocmanager r-ggplot2 r-seurat=3.2.3 cmake

