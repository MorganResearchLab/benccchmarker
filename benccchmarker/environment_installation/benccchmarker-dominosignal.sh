mamba create -n benccchmarker-dominosignal -c conda-forge -c bioconda -c defaults python=3.8 r-base=4.3.3 r-devtools r-biocmanager r-ggplot2 r-seurat cmake scanpy loompy anndata numpy -y

mamba activate benccchmarker-dominosignal

echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

Rscript -e "options(repos = c(PkgMgr='https://packagemanager.rstudio.com/all/__linux__/focal/latest'))"

Rscript -e "devtools::install_github('FertigLab/dominoSignal')"

# Python library requirements for generating PyScenic input
pip install pyscenic

wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/refs/heads/master/example/allTFs_hg38.txt
wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/refs/heads/master/example/motifs.tbl

wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather


