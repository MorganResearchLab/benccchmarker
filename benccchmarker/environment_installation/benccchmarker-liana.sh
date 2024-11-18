mamba create -n benccchmarker-liana -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.2.0 -y

mamba activate benccchmarker-liana

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

pip install -y liana