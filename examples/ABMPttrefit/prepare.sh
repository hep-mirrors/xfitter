#!/bin/bash

if [ ! -d datafiles/fixedTarget/nomad ]; then
rm datafiles
wget "https://cernbox.cern.ch/remote.php/dav/public-files/grDRAZ7VnuwIjYw/abmp16-newdata.tar.gz"
tar xvzpf abmp16-newdata.tar.gz
rm abmp16-newdata.tar.gz
fi

if [ ! -d datafiles/lhc/cms/topProduction/2303.15340 ]; then
wget "https://cernbox.cern.ch/remote.php/dav/public-files/aGXPdhzBXa9yH4y/abmptt-newdata.tar.gz"
tar xvzpf abmptt-newdata.tar.gz
rm abmptt-newdata.tar.gz
fi

if [ ! -d ../../fewzgrids ]; then
wget "https://cernbox.cern.ch/remote.php/dav/public-files/CnmQytk8SYDCMJi/fewzgrids.tar.gz"
tar xvzpf fewzgrids.tar.gz -C ../..
rm fewzgrids.tar.gz
ln -s `pwd`/../../fewzgrids .
else
ln -s ../../fewzgrids
fi
