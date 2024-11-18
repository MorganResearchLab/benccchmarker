mamba create -n benccchmarker-icellnetv1 -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 -y

mamba activate benccchmarker-icellnetv1

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install -y r-devtools r-biocmanager r-ggplot2 r-seurat cmake

Rscript -e 'install.packages(c("devtools", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist"))'
Rscript -e 'BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate", "jetset", "AnnotationDbi", "hgu133plus2.db"))'
Rscript -e "devtools::install_github('soumelis-lab/ICELLNET', ref = 'v1.3.0', subdir='icellnet')"