# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/JiangBioLab/DeepCCI.git

mamba create -n benccchmarker-deepcci -c conda-forge -c bioconda -c defaults python=3.7.4 r-base pip -y

mamba activate benccchmarker-deepcci

cd DeepCCI

mamba install -y rust

pip install -r requirements.txt

Rscript -e install.packages('Seurat', dependencies=TRUE)
Rscript -e install.packages('igraph', dependencies=TRUE)
Rscript -e install.packages('NMF', dependencies=TRUE)
Rscript -e install.packages('devtools', dependencies=TRUE)
Rscript -e devtools::install_github('jokergoo/circlize')
Rscript -e devtools::install_github('jokergoo/ComplexHeatmap')
Rscript -e devtools::install_github('sqjin/CellChat')
Rscript -e devtools::install_github('satijalab/seurat-data')
