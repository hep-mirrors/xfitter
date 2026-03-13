
-- A simple fit of NLO pion PDFs based on the Fantomas 1.0 parametrization
    https://arxiv.org/abs/2505.13594


### How to run ###
git clone  ssh://git@gitlab.cern.ch:7999/fitters/xfitter.git

cd xfitter
git checkout Fantomas
./make.sh install
cd examples/Pion_Fantomas/
ln -s ~/YourDataFile/datafiles .
../../bin/xfitter 
