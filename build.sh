#!/bin/bash

# unload potentially conflicting anaconda instances
{ # try
    module unload python-anaconda2 &&
    module unload python-anaconda3
} || { # catch
    echo 'module unloading failed: maybe module does not exist'
}


# install miniconda for local independent package management
wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh -O miniconda.sh
mkdir dependencies
chmod 775 miniconda.sh
chmod 775 dependencies
bash miniconda.sh -b -p ./dependencies/miniconda
rm miniconda.sh

# activate conda virtual environment
source ./dependencies/miniconda/bin/activate

conda install click

# install pyflow for automated task management
pip install https://github.com/Illumina/pyflow/releases/download/v1.1.17/pyflow-1.1.17.tar.gz

conda create -n py3 python=3.6

source activate py3
pip install requirements.txt
conda install -c r r
conda install click
Rscript install.r 
