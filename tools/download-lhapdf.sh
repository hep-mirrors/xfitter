#!/bin/bash

# manual script to download LHAPDF sets needed by some xfitter tests, 
# it is used by .gitlab-ci.yml, because `lhapdf install` seems to be very tricky

LINK='http://lhapdfsets.web.cern.ch/lhapdfsets/current/'

for pdf in $@; do
  wget ${LINK}'/'${pdf}'.tar.gz'
  olddir=`pwd`
  cd `lhapdf-config --datadir`
  tar xvzpf ${olddir}'/'${pdf}'.tar.gz'
  cd - > /dev/null
  rm ${pdf}'.tar.gz'
done
  
  
