# Requires its own environment

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src

git clone https://github.com/QSong-github/spaCI.git
cd spaCI

mamba create -n benccchmarker-spaci -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 -y 

mamba activate benccchmarker-spaci

pip install git+https://github.com/cailab-tamu/spaci.git 
