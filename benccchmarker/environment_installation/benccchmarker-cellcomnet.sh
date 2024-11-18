cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src

git clone https://github.com/ringsssss/CellComNet.git

mamba create -n benccchmarker-cellcomnet -c conda-forge -c bioconda -c defaults python=3.10 -y

mamba activate benccchmarker-cellcomnet

mamba install -y tensorflow=1.14.0 pandas=1.1.5 numpy scikit-learn

pip install keras anndata scanpy