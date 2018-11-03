#!/usr/bin/env bash

CONDAENV=phenoX

if ! (which wget); then
    echo "No `wget` installation found. Please install wget first."
    exit 1
fi

if ! (which conda); then
	echo "No `conda` installation found.  Installing..."
	if [[ $(uname) == "Darwin" ]]; then
	  wget --continue http://repo.continuum.io/archive/Anaconda3-4.4.0-MacOSX-x86_64.sh
	  bash Anaconda3-4.4.0-MacOSX-x86_64.sh -b
	else
	  wget --continue http://repo.continuum.io/archive/Anaconda3-4.4.0-Linux-x86_64.sh
	  bash Anaconda3-4.4.0-Linux-x86_64.sh -b
	fi
    printf '\n# Anaconda path addition\n' >> $HOME/.bashrc
    echo 'export PATH="$HOME/anaconda3/bin:$PATH"' >> $HOME/.bashrc
fi

CONDAPATH="$(which conda)" || CONDAPATH=$HOME/anaconda3/bin
echo "CONDAPATH => $CONDAPATH"

export PATH=$CONDAPATH:$PATH

conda create -n ${CONDAENV} -y python==3.6 pip pytest || true

echo "Activating Conda Environment ----->"
conda activate ${CONDAENV} || source $CONDAPATH/activate ${CONDAENV} 

echo "Installing required Python libraries ----->"
pip install -r requirements.txt
conda install -y matplotlib

echo "Installing required conda-forge libraries ----->"
conda install -y -c conda-forge rpy2 r-ape r-pvclust r-circlize readline r-gplots

echo "Installing phenoX ----->"
python setup.py develop

echo "Installing SpaCy en ----->"
python -m spacy download en
