# Requires its own environment

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src

git clone https://github.com/cailab-tamu/scTenifoldXct.git
cd scTenifoldXct
sed -i 's/name: scTenifold/name: benccchmarker-sctenifoldxct/' /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src/scTenifoldXct/environment.yml
mamba env create -f environment.yml

mamba activate benccchmarker-sctenifoldxct

pip install git+https://github.com/cailab-tamu/scTenifoldXct.git 
