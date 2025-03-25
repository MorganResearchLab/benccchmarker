# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/SysBioOncology/RaCInG.git

mamba create -n benccchmarker-racing -c conda-forge -c bioconda -c defaults python=3.8.11 r-base=4.2.0 -y

mamba activate benccchmarker-racing

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install bioconductor-omnipathr r-immunedeconv bioconductor-easier r-epic r-mcpcounter bioconductor-quantiseqr r-xcell r-consensustme r-corrplot r-dplyr r-ggplot2 -y

pip install liana numpy==1.16.6 scipy==1.7.1 matplotlib==3.4.3 seaborn==0.11.2 pandas==1.2.4 Adjusttext==0.7.3