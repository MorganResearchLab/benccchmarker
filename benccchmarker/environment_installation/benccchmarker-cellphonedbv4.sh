mamba create -n benccchmarker-cellphonedbv4 -c conda-forge -c bioconda -c defaults python=3.10 -y

mamba activate benccchmarker-cellphonedbv4

pip install git+https://github.com/ventolab/CellphoneDB.git@v4.0.0