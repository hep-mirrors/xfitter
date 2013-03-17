#!/bin/sh
#
echo "PBS_O_HOST=$PBS_O_HOST"
echo "PBS_O_HOME=$PBS_O_HOME"
echo "PBS_O_WORKDIR=$PBS_O_WORKDIR"
module load gcc
PROGRAM=FitPDF
#
LOGFILE=${PROGRAM}.log
ERRFILE=${PROGRAM}.err
PROGRAM_EXE=${PROGRAM}

cd $PBS_O_WORKDIR

if [ -n "$RunTimeDir" ]
then
  RTdir=$RunTimeDir/$PBS_JOBID
  mkdir -p $RTdir
  cp -r * $RTdir
  cd $RTdir
fi

unzip input.zip
rm input.zip
mkdir -p output
#
echo  "Running program $PROGRAM ........................."
#
time ./$PROGRAM_EXE > $LOGFILE 2> $ERRFILE
#
#  --- Check status
if [ $? != 0 ]
then
  echo "ERROR running $PROGRAM."
  exit 1
fi
echo "Successfully ran $PROGRAM."
echo ''
ls -l
echo ''
echo ' ------------------------------------------------------- '
echo ''
# ls -l $PBS_O_WORKDIR
# echo ''
# echo ' ------------------------------------------------------- '
# echo ''
rm -f *.wgt
rm -f $PROGRAM_EXE
zip -jo output.zip output/*
if [ -n "$RunTimeDir" ]
then
  echo ''
  echo "Copying files to $PBS_O_WORKDIR ........"
  echo ''
  cp -p * $PBS_O_WORKDIR
fi
date
echo "$PROGRAM finished."
