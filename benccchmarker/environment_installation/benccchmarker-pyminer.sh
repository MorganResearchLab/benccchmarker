mamba create -n benccchmarker-pyminer -c conda-forge -c bioconda -c defaults python=3.8 -y

mamba activate benccchmarker-pyminer

pip install numpy
mamba install -y "setuptools <65"
pip install bio-pyminer