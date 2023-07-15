#!/bin/bash

# manual script to download LHAPDF sets needed by some xfitter tests, 
# it is used by .gitlab-ci.yml, because `lhapdf install` seems to be very tricky

LINK='http://lhapdfsets.web.cern.ch/lhapdfsets/current/'

# check if wget available, otherwise revert to curl
 which wget >& /dev/null
if [[ $? == 0 ]]
then
    http=wget
else
    which curl >& /dev/null
    if [[ $? == 0 ]]
    then
        http=curl
    else
        echo "Error, wget or curl not found"
        exit
    fi
fi


# download and move to datadir list of input PDFs
for pdf in $@; do
  if [[ $http == "curl" ]]
  then
    curl -L -O ${LINK}'/'${pdf}'.tar.gz'
  else 
    wget ${LINK}'/'${pdf}'.tar.gz'
  fi
  olddir=`pwd`
  cd `lhapdf-config --datadir`
  tar xvzpf ${olddir}'/'${pdf}'.tar.gz'
  cd - > /dev/null
  rm ${pdf}'.tar.gz'
done
  
  
