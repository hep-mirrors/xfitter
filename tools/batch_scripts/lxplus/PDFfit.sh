#!/bin/bash

outName="epWZ16+ttbar"                #Output Name
emailName="francesco.giuli@cern.ch"   #User email

# send email to know when the job started
echo "PDF ${outName} fit has started" | mail -s "PDF fit batch/queue status" ${emailName}

# creates a local log every time the script runs
cd /afs/cern.ch/work/f/fgiuli ; echo " " >> submit_hist.txt ; echo "Fit on ${outName}" >> submit_hist.txt

# load all the necessary PAHT, LD_LIBRARY_FLAGS, CXXFLAGS, etc.
cd /afs/cern.ch/work/f/fgiuli/public ; source .bashrc;
cd /afs/cern.ch/work/f/fgiuli/public/installation/xfitter ; source setup.sh

# check that files were sourced correctly, write it in the local log
cd /afs/cern.ch/work/f/fgiuli/public ; echo ">>> .bashrc sourced, setup.sh sourced, running xfitter located at:" >> submit_hist.txt ; which xfitter >> submit_hist.txt ; date >> submit_hist.txt

# location of the steering file
cd /afs/cern.ch/work/f/fgiuli/public/globalFit

xfitter

# send second email, whether the fit was completed succesfully or failed
if [[ $? != 0 ]]
then
    cd /afs/cern.ch/work/f/fgiuli/public ; echo ">>>>>> fit finished with an error:" >> submit_hist.txt ; date >> submit_hist.txt ; echo " " >> submit_hist.txt
    echo "ERROR: PDF ${outName} fit stopped unexpectedly" | mail -s "PDF fit batch/queue status" ${emailName}
    exit
fi
echo "PDF ${outName} fit has finished" | mail -s "PDF fit batch/queue status" ${emailName}

cd /afs/cern.ch/work/f/fgiuli/public ; echo ">>> fit finished:" >> submit_hist.txt ; date >> submit_hist.txt ; echo " " >> submit_hist.txt
