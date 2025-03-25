mamba create -n benccchmarker-scseqcomm -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 r-devtools r-biocmanager r-ggplot2 r-seurat cmake r-magick bioconductor-clusterprofiler r-peakram -y

mamba activate benccchmarker-scseqcomm

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

Rscript -e "install.packages('add2ggplot', dependencies = TRUE)"
Rscript -e 'devtools::install_gitlab("sysbiobig/scseqcomm")'

