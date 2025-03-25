mamba create -n benccchmarker-pyminer -c conda-forge -c bioconda -c defaults python=3.8 pillow networkx seaborn numpy scikit-learn h5py matplotlib scipy gprofiler-official umap-learn python-louvain statsmodels -y

mamba activate benccchmarker-pyminer

pip install numpy scanpy
mamba install -y "setuptools <65"
pip install bio-pyminer

#rsync -av --exclude='__init__.py' --exclude='__pycache__' /uoa/scratch/users/r04mr23/.conda/envs/benccchmarker-pyminer/lib/python3.8/site-packages/pyminer_norm/*.py /uoa/scratch/users/r04mr23/.conda/envs/benccchmarker-pyminer/lib/python3.8/site-packages/pyminer/