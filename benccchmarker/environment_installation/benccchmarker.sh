mamba create -n benccchmarker -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.1

mamba activate benccchmarker
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

mamba install -y r-devtools r-biocmanager r-ggplot2 r-seurat cmake rpy2 r-peakram r-anndata

# Install CellChat
Rscript -e "BiocManager::install('Biobase', dependencies = TRUE)"
mamba install -y r-ggplot2 r-nmf r-presto r-rspectra
Rscript -e "BiocManager::install('ComplexHeatmap', dependencies = TRUE)"
Rscript -e "BiocManager::install('BiocNeighbors', dependencies = TRUE)"
Rscript -e "install.packages('ggpubr', dependencies = TRUE)"
Rscript -e "devtools::install_github('jinworks/CellChat')"

# Install CellCall
Rscript -e "BiocManager::install('DOSE', dependencies = TRUE)"
Rscript -e "BiocManager::install('enrichplot', dependencies = TRUE)"
Rscript -e "BiocManager::install('clusterProfiler', dependencies = TRUE)"
Rscript -e "devtools::install_github('ShellyCoder/cellcall')"

# Install celltalker
Rscript -e "devtools::install_github('arc85/celltalker')"

# Install crosstalker
Rscript -e "devtools::install_github('https://github.com/CostaLab/CrossTalkeR')"

# Install iCellNet v2
Rscript -e 'install.packages(c("devtools", "jetset", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer", "data.table", "grid", "gridExtra", "ggthemes", "scales","rlist"))'
Rscript -e 'BiocManager::install(c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate"))'
Rscript -e 'devtools::install_github("soumelis-lab/ICELLNET",ref="master", subdir="icellnet")'

# Install scMLnet
Rscript -e "install.packages('parallel', dependencies = TRUE)"
Rscript -e 'devtools::install_github("YUZIXD/scMLnet")'

# Install SingleCellSignalR
Rscript -e "BiocManager::install('SingleCellSignalR', dependencies = TRUE)"

# Install RSoptSC
Rscript -e 'devtools::install_github("mkarikom/RSoptSC")'

# Install scSeqComm
mamba install -y r-magick
Rscript -e "install.packages('add2ggplot', dependencies = TRUE)"
Rscript -e 'devtools::install_gitlab("sysbiobig/scseqcomm")'

# Install scDiffCom
Rscript -e 'devtools::install_github("CyrilLagger/scDiffCom")'

# Install NicheNet
Rscript -e 'devtools::install_github("saeyslab/nichenetr")'

# Install multinichenetr
Rscript -e 'devtools::install_github("saeyslab/multinichenetr")'

# Install CytoTalk
pip install git+https://github.com/fraenkel-lab/pcst_fast.git
Rscript -e 'devtools::install_github("tanlabcode/CytoTalk")'

# Install cellphonedb
pip install cellphonedb==5.0.0





## FAILS NEED NEW ENVIRONMENT

# Install iTALK (Need Cairo)
Rscript -e "BiocManager::install('monocle', dependencies = TRUE)"
mamba install r-cairo cairo
Rscript -e "BiocManager::install('scater', dependencies = TRUE)"
## Rscript -e "install.packages('Cairo', dependencies = TRUE)"
Rscript -e "BiocManager::install('scde', dependencies = TRUE)"
Rscript -e "devtools::install_github('Coolgenome/iTALK', build_vignettes = TRUE)"

# Install COMUNET
mamba install r-sdmtools
Rscript -e "devtools::install_github('ScialdoneLab/COMUNET/COMUNET')"