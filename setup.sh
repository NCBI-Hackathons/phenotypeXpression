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

export PATH=$HOME/anaconda3/bin:$PATH

conda create -n ${CONDAENV} -y python==3.6 pip pytest || true

echo "Activating Conda Environment ----->"
source ~/anaconda3/bin/activate ${CONDAENV}

pip install -r requirements.txt

python setup.py develop