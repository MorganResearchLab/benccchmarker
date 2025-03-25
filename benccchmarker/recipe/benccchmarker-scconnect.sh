mamba create -n benccchmarker-scconnect -c conda-forge -c bioconda -c defaults python=3.9 pip -y

mamba activate benccchmarker-scconnect

pip install git+https://github.com/JonETJakobsson/scConnect

# Change iteritems to items on connect.py
pip install pandas==1.5.3 
