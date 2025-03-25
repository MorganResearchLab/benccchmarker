cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/yuanruya/CellEnBoost.git

mamba create -n benccchmarker-cellenboost -c conda-forge -c bioconda -c defaults python=3.10 pip -y

mamba activate benccchmarker-cellenboost

pip install tensorflow==1.14.0 keras pandas numpy==1.24.3 scikit-learn lightgbm

cd cellenboost