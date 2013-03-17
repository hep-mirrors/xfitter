#!/bin/sh
#
#   Turn on Verbose Mode
set -x
#
PROGRAM=FitPDF
#
LOGFILE=${PROGRAM}.log
ERRFILE=${PROGRAM}.err
PROGRAM_EXE=${PROGRAM}

unzip input.zip
rm input.zip
mkdir -p output
#
echo  "Running program $PROGRAM ........................."
#
time ./$PROGRAM_EXE > $LOGFILE 2> $ERRFILE
#
#   Check the status
#   ----------------
if [ $? != 0 ]
then
        echo "ERROR running $PROGRAM."
        exit 1
fi
echo "Successfully ran $PROGRAM."
echo ''
#
ls -l
#
echo ''
echo ' ------------------------------------------------------- '
echo ''
zip -jo output.zip output/*
date
echo "$PROGRAM Done."
#
#   Turn off Verbose Mode
set -
#
