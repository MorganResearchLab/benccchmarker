cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src
git clone https://github.com/Xwhut/CellDialog.git

mamba create -n benccchmarker-celldialog -c conda-forge -c bioconda -c defaults python=3.9 pip -y

mamba activate benccchmarker-celldialog

pip install seaborn>=0.12.0 pandas>=1.5.3 numpy>=1.25.2 scikit-learn>=1.2.1 matplotlib>=3.7.2 xgboost==1.6.2 lightgbm==3.3.0 KTBoost==0.2.2 gpboost==0.7.10 anndata


cd CellDialog