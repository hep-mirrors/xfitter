#!/bin/sh
#
#   Turn on Verbose Mode
#set -x

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

# This is SUFFIX for input file names in $INPUTDIR: could be e.g. 'def' (by default), 'dipole', etc.
SUFFIX=$1
SUFFIXDEF='def'
if [ -z $SUFFIX ]; then
  SUFFIX=$SUFFIXDEF
fi

rm -rf temp
mkdir temp

echo "========================================"
echo "Running checks"
echo "========================================"
echo "validation test: $SUFFIX"
echo "PASS if code runs properly"
echo "FAIL if code fails to reproduce expected results"
echo "========================================"

INPUTDIR="input_steering_22"
EXAMPLEDIR="examples_22/output-$SUFFIX"

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
else
echo "========================================"
echo "Failed validation with default steering"
echo "========================================"
flagAllFine=0
fi 
echo "Checking all output files ..."
for file in `find output -type f`; do
  targetfile=${EXAMPLEDIR}/`echo ${file} | sed -e 's/output//'`
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
if [ $flagAllFine != 0 ]; then
  echo "Everything is fine"
else
  echo "Something failed: see above for details"
fi
echo "========================================"

rm -rf temp

if [ $flagAllFine = 0 ]; then
  exit 1
fi

