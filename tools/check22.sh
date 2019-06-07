#!/bin/sh
#
#   Turn on Verbose Mode
#set -x

# function to check if $1 and $2 are identical
checkFile()
{
  printf "diff $1 $2 ... "
  diff $1 $2 > /dev/null
  exitcode=$?

  if [ $exitcode = 0 ]; then
  echo "OK"
  else
  echo "FAILED"
  fi 
}

# the scrpit starts here
if [ "$#" -eq 0 ] || [ "$1" = "--help" ]; then
  echo "Usage: check.sh <TEST> [OPTION]"
  echo "OPTION:"
  echo "  --copy to copy results and make them reference"
  echo "  --help to see this message"
  exit 1
fi

# This is SUFFIX for input file names in $INPUTDIR: could be e.g. 'def' (by default), 'dipole', etc.
SUFFIX=$1
SUFFIXDEF='def'
if [ -z $SUFFIX ]; then
  SUFFIX=$SUFFIXDEF
fi

# Optionally copy results and make them reference
COPYRESULTS=0
if [ "$2" = "--copy" ]; then
  COPYRESULTS=1
fi

rm -rf temp
mkdir temp

echo "========================================"
echo "Running checks: $SUFFIX"
echo "========================================"
if [ $COPYRESULTS -eq 0 ]; then
  echo "validation test:"
  echo "PASS if code runs properly"
  echo "FAIL if code fails to reproduce expected results"
  echo "========================================"
fi

INPUTDIR="input_steering_22"
EXAMPLEDIR="examples_22/output-$SUFFIX"
if [ $COPYRESULTS -eq 1 ]; then
  rm -rf $EXAMPLEDIR
  mkdir -p $EXAMPLEDIR
  mkdir -p $EXAMPLEDIR'/xfitter_pdf'
  echo "Results will be stored as reference in $EXAMPLEDIR"
  echo "========================================"
fi

FlagUnique=0
STEERING=${INPUTDIR}/steering.txt.${SUFFIX}
if [ ! -f $STEERING ]; then
  STEERING=${INPUTDIR}/steering.txt.${SUFFIXDEF}
else
  FlagUnique=1
fi
cp ${STEERING} steering.txt

PARAMETERS=${INPUTDIR}/parameters.yaml.${SUFFIX}
if [ ! -f $PARAMETERS ]; then
  PARAMETERS=${INPUTDIR}/parameters.yaml.${SUFFIXDEF}
else
  FlagUnique=1
fi
cp ${PARAMETERS} parameters.yaml

CONSTANTS=${INPUTDIR}/constants.yaml.${SUFFIX}
if [ ! -f $CONSTANTS ]; then
  CONSTANTS=${INPUTDIR}/constants.yaml.${SUFFIXDEF}
else
  FlagUnique=1
fi
cp ${CONSTANTS} constants.yaml

if [ $FlagUnique = 0 ] && [ $SUFFIX != "def" ]; then
  echo "Failed to find input files for test \"$SUFFIX\""
  exit 1
fi

echo "Using ${STEERING}"
echo "Using ${PARAMETERS}"
echo "Using ${CONSTANTS}"
echo "========================================"

bin/xfitter >/dev/null

flagAllFine=0
if [ $COPYRESULTS -eq 0 ]; then
  grep  'After' output/Results.txt > temp/out.txt
  grep  'After' ${EXAMPLEDIR}/Results.txt > temp/def.txt

  cat temp/out.txt
  diff temp/out.txt temp/def.txt 

  exitcode=$?

  flagAllFine=1
  if [ $exitcode = 0 ]; then
    echo "========================================"
    echo "Check of chi^2 is fine"
    echo "========================================"
    flagAllFine=1
  else
    echo "========================================"
    echo "Failed validation with default steering"
    echo "========================================"
    flagAllFine=0
  fi
fi

echo "Checking all output files ..."
for file in `find output -type f`; do
  targetfile=${EXAMPLEDIR}/`echo ${file} | sed -e 's/output//'`
  if [ $COPYRESULTS -eq 1 ]; then
    echo "storing $file as $targetfile"
    cp $file $targetfile
    continue
  fi
  if [ ! -f $targetfile ]; then
    echo "No reference $targetfile"
    flagAllFine=0
    continue
  fi
  out=`checkFile $file $targetfile`
  echo $out
  echo $out | grep FAILED > /dev/null
  if [ $exitcode = 1 ]; then
    flagAllFine=0
  fi
done
echo "========================================"
if [ $COPYRESULTS -eq 0 ]; then
  if [ $flagAllFine != 0 ]; then
    echo "Everything is fine"
  else
    echo "Something failed: see above for details"
  fi
  echo "========================================"
fi

rm -rf temp

if [ $flagAllFine = 0 ]; then
  exit 1
fi
