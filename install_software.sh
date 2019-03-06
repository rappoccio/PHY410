# Ubuntu+tools, python, pip, python modules
apt-get update
apt-get install -y wget g++ libtool rsync make x11-apps python-dev python-numpy python-pip python-tk 
rm -rf /var/lib/apt/lists/*
#pip install --no-cache-dir matplotlib scipy numpy scikit-learn keras tensorflow jupyter metakernel zmq notebook==5.* plaidml-keras plaidbench energyflow
pip2 install --upgrade pip && pip2 install --no-cache-dir matplotlib==2.1 scipy numpy scikit-learn keras tensorflow ipykernel==4.10.0 ipython==5.7 jupyter metakernel zmq notebook==5.*
make



