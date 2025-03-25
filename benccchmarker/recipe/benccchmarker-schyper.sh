mamba create -n benccchmarker-schyper -c conda-forge -c bioconda -c defaults python=3.10 r-base=4.3.3 r-devtools -y

mamba activate benccchmarker-schyper

git clone https://github.com/Lwyonly/scHyper.git

pip install -e git+https://github.com/Lwyonly/scHyper.git

mamba install pytorch==2.5.1 torchvision==0.20.1 torchaudio==2.5.1  pytorch-cuda=11.8 -c pytorch -c nvidia -y

# No such file or directory: '/uoa/scratch/users/r04mr23/.conda/envs/benccchmarker-schyper/lib/python3.10/site-packages/scHyper/database/human.csv'

# (benccchmarker-schyper) [r04mr23@maxlogin2(maxwell) benccchmarker-schyper]$ mkdir /uoa/scratch/users/r04mr23/.conda/envs/benccchmarker-schyper/lib/python3.10/site-packages/scHyper/database/
# (benccchmarker-schyper) [r04mr23@maxlogin2(maxwell) benccchmarker-schyper]$ cp /uoa/home/r04mr23/sharedscratch/trial/benCCChmarker/benccchmarker/trial/benccchmarker-schyper/scHyper/scHyper/database/human.csv /uoa/scratch/users/r04mr23/.conda/envs/benccchmarker-schyper/lib/python3.10/site-packages/scHyper/database/human.csv
# (benccchmarker-schyper) [r04mr23@maxlogin2(maxwell) benccchmarker-schyper]$ 