# Requires Seurat v3 and some packages that will clash with the default benccchmarker-all packages

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/forrest-lab/NATMI.git

mamba create -n benccchmarker-natmi -c conda-forge -c bioconda -c defaults python=3.7.6 r-base pip pandas xlsxwriter xlrd igraph seaborn networkx bokeh holoviews openpyxl anndata scanpy pygraphviz -y

mamba activate benccchmarker-natmi

pip install "h5py<3.2"

cd natmi