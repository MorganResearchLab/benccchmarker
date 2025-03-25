# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

mamba create -n benccchmarker-rnamagnet -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 -y

mamba activate benccchmarker-rnamagnet

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install -y r-devtools r-biocmanager r-ggplot2 r-seurat=3.2.3 cmake

# Install RNAMagnet
Rscript -e 'devtools::install_github("veltenlab/rnamagnet")'