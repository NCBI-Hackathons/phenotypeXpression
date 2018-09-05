CONDAENV=phenoX

if ! (which wget); then
    echo "No `wget` installation found. Please install wget first."
    exit 1
fi

if ! (which conda); then
	echo "No `conda` installation found.  Installing..."
	if [[ $(uname) == "Darwin" ]]; then
	  wget --continue http://repo.continuum.io/archive/Anaconda3-4.3.1-MacOSX-x86_64.sh
	  bash Anaconda3-4.3.1-MacOSX-x86_64.sh -b
	else
	  wget --continue http://repo.continuum.io/archive/Anaconda3-4.3.1-Linux-x86_64.sh
	  bash Anaconda3-4.3.1-Linux-x86_64.sh -b
	fi
fi

CONDAPATH="$(which conda)"

export PATH=$CONDAPATH:$PATH

conda create -n ${CONDAENV} -y python==3.6 pip pytest || true

echo "Activating Conda Environment ----->"
conda $CONDAPATH/activate ${CONDAENV}

echo "Installing required Python libraries ----->"
pip install -r requirements.txt
conda install matplotlib

echo "Installing required R libraries ----->"
conda install -c r r-base r-ape
conda install rpy2
conda install -c conda-forge r-pvclust r-circlize

echo "Installing phenoX ----->"
python setup.py develop

echo "Installing SpaCy en ----->"
python -m spacy download en
