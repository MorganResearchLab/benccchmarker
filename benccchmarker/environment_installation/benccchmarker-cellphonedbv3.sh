mamba create -n benccchmarker-cellphonedbv3 -c conda-forge -c bioconda -c defaults python=3.10 -y

mamba activate benccchmarker-cellphonedbv3

pip install git+https://github.com/ventolab/CellphoneDB.git@v3.0.0