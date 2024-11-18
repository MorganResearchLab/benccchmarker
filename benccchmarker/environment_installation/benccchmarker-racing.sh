# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/SysBioOncology/RaCInG.git

mamba create -n benccchmarker-racing -c conda-forge -c bioconda -c defaults python=3.8.11 r-base=4.2.0 -y

mamba activate benccchmarker-racing

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

pip install liana

mamba install bioconductor-omnipathr r-immunedeconv bioconductor-easier epic r-mcpcounter bioconductor-quantiseqr r-xcell r-consensustme r-corrplot r-dplyr r-ggplot2