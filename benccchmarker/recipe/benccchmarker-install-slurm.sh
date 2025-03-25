#!/bin/bash
## SBATCH --time=03:00:00

#SBATCH --ntasks=1
#SBATCH --mem=8g
#SBATCH --job-name=benccchmarker-environment-install
#SBATCH --output=benccchmarker-environment-install.log

module load mamba
source ~/.bash_profile
source ~/.bashrc

mamba activate benccchmarker

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