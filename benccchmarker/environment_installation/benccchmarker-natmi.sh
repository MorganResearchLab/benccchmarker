# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/forrest-lab/NATMI.git

mamba create -n benccchmarker-natmi -c conda-forge -c bioconda -c defaults python=3.7.6 r-base pip -y

mamba activate benccchmarker-natmi
mamba install pygraphviz -y

pip install pandas==1.0.3 XlsxWriter==1.2.8 xlrd==1.2.0 igraph seaborn==0.10.1 igraph NetworkX==2.4 bokeh==2.0.2 and holoviews==1.13.2

cd natmi