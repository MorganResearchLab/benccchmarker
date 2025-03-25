# Requires its own environment

cd /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src

git clone https://github.com/miladrafiee/disir_package.git
cd disir_package
sed -i 's/name: disir_package/name: benccchmarker-disir/' /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/src/disir_package/environment.yml
mamba env create -f environment.yml

mamba activate benccchmarker-disir

pip install .
