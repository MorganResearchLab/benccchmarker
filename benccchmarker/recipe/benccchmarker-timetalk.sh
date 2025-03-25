# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

mamba create -n benccchmarker-timetalk -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 -y

mamba activate benccchmarker-timetalk

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install -y r-devtools r-biocmanager r-ggplot2 r-seurat cmake

mamba install -y r-monocle3
Rscript -e 'devtools::install_github("shenorrLab/cellAlign")'
Rscript -e 'BiocManager::install("RTN", dependencies = TRUE)'
Rscript -e 'devtools::install_github("ChengLiLab/TimeTalk",ref="main")'