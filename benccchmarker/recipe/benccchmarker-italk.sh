mamba create -n benccchmarker-italk -c conda-forge -c bioconda -c defaults r-base=4.3.3 r-devtools r-biocmanager r-ggplot2 r-seurat cmake r-cairo cairo r-peakram -y

mamba activate benccchmarker-italk

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

Rscript -e "BiocManager::install('monocle', dependencies = TRUE)"
Rscript -e "BiocManager::install('scater', dependencies = TRUE)"
Rscript -e "BiocManager::install('scde', dependencies = TRUE)"
Rscript -e "devtools::install_github('Coolgenome/iTALK', build_vignettes = TRUE)"