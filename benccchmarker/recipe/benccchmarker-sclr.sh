mamba create -n benccchmarker-sclr -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 -y

mamba activate benccchmarker-sclr

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install -y r-devtools r-biocmanager r-ggplot2 r-seurat cmake

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/cyhsuTN/scLR.git

cd scLR