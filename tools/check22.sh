#!/bin/bash
#
#   Turn on Verbose Mode
#set -x

# output style
FAILED="\e[31m\e[1mFAILED\e[0m" # bold red
OK="\e[32m\e[1mOK\e[0m" # bold green

# function to check if $1 and $2 are identical
checkFile()
{
  printf "diff $1 $2 ... "
  diff $1 $2 > /dev/null
  exitcode=$?

  if [ $exitcode = 0 ]; then
  echo -e "$OK"
  else
  echo -e "$FAILED"
  fi 
}

# the scrpit starts here
if [ "$#" -eq 0 ] || [ "$1" = "--help" ]; then
  echo "Usage: check.sh <TEST> [OPTION]"
  echo "OPTION could be:"
  echo "  --copy to copy results and make them reference"
  echo "  --help to see this message"
  exit 1
fi

# TESTNAME is the name of directory in examples_22/
TESTNAME=$1
TESTNAMEDEF='defaultNNLO'
if [ -z $TESTNAME ]; then
  TESTNAME=$TESTNAMEDEF
fi

# Optionally copy results and make them reference
COPYRESULTS=0
if [ "$2" = "--copy" ]; then
  COPYRESULTS=1
fi

rm -rf temp/$TESTNAME
mkdir -p temp/$TESTNAME

echo "========================================"
echo "Running checks: $TESTNAME"
echo "========================================"
if [ $COPYRESULTS -eq 0 ]; then
  echo "validation test:"
  echo "OK if code runs properly"
  echo "FAILED if code fails to reproduce expected results"
  echo "========================================"
fi

INPUTDIR="examples_22/$TESTNAME"
EXAMPLEDIR="examples_22/$TESTNAME/output"
if [ $COPYRESULTS -eq 1 ]; then
  rm -rf $EXAMPLEDIR
  mkdir -p $EXAMPLEDIR
  #mkdir -p $EXAMPLEDIR'/xfitter_pdf'
  echo "Results will be stored as reference in $EXAMPLEDIR"
  echo "========================================"
else
  echo "Output will be stored in temp/$TESTNAME/xfitter.txt"
  echo "========================================"
fi

if [ ! -d $INPUTDIR ]; then
  echo "Failed to find input files for test \"$TESTNAME\""
  echo "expected directory $INPUTDIR"
  exit 1
fi
echo "Using input files from ${INPUTDIR}"
echo "========================================"

cp ${INPUTDIR}/steering.txt ./
cp ${INPUTDIR}/parameters.yaml ./
cp ${INPUTDIR}/constants.yaml ./

bin/xfitter >& temp/$TESTNAME/xfitter.txt

flagBAD=0
if [ $COPYRESULTS -eq 0 ]; then
  grep  'After' output/Results.txt > temp/out.txt
  grep  'After' ${EXAMPLEDIR}/Results.txt > temp/def.txt

  cat temp/out.txt
  diff temp/out.txt temp/def.txt 

  exitcode=$?
  rm -f temp/out.txt temp/def.txt

  if [ $exitcode = 0 ]; then
    echo "========================================"
    echo -e "Check of chi^2 is $OK"
    echo "========================================"
  else
    echo "========================================"
    echo -e "$FAILED validation with default steering"
    echo "========================================"
    flagBAD=1
  fi
fi

echo "Checking all output files ..."
for file in `find output -type f`; do
  targetfile=${EXAMPLEDIR}/`echo ${file} | sed -e 's/output//'`
  if [ $COPYRESULTS -eq 1 ]; then
    echo "storing $file as $targetfile"
    if [ ! -d `dirname $targetfile` ]; then
      mkdir `dirname $targetfile`
    fi
    cp $file $targetfile
    continue
  fi
  if [ ! -f $targetfile ]; then
    echo "No reference $targetfile"
    flagBAD=1
    continue
  fi
  out=`checkFile $file $targetfile`
  echo $out
  echo $out | grep FAILED > /dev/null
  exitcode=$?
  if [ $exitcode == 0 ]; then
    flagBAD=1
  fi
done
echo "========================================"
if [ $COPYRESULTS -eq 0 ]; then
  if [ $flagBAD == 0 ]; then
    echo -e "Everything is $OK"
  else
    echo -e "Something $FAILED: see above for details"
  fi
  echo "========================================"
fi

exit $flagBAD
