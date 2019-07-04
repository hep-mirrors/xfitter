#!/bin/bash
#
#   Turn on Verbose Mode
#set -x

# output style
FAILED="\e[31m\e[1mFAILED\e[0m" # bold red
OK="\e[32m\e[1mOK\e[0m" # bold green

# xfitter binary
xfitter=`pwd`/'bin/xfitter'

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

rundir=temp/$TESTNAME
rm -rf $rundir
mkdir -p $rundir

echo "========================================"
echo "Running check: $TESTNAME"
echo "========================================"
if [ $COPYRESULTS -eq 0 ]; then
  echo "This is validation test:"
  echo "OK if code runs properly"
  echo "FAILED if code fails to reproduce expected results"
  echo "========================================"
fi

INPUTDIR="examples_22/$TESTNAME"
EXAMPLEDIR="examples_22/$TESTNAME/output"
if [ $COPYRESULTS -eq 1 ]; then
  echo "Results will be stored as reference in $EXAMPLEDIR"
fi
echo "Running in temp/$TESTNAME"
echo "Log file stored in temp/$TESTNAME/xfitter.txt"

if [ ! -d $INPUTDIR ]; then
  echo "Failed to find input files for test \"$TESTNAME\""
  echo "expected directory $INPUTDIR"
  exit 1
fi
echo "Using input files from ${INPUTDIR}"
echo "========================================"

cp ${INPUTDIR}/steering.txt $rundir
cp ${INPUTDIR}/parameters.yaml $rundir
cp ${INPUTDIR}/constants.yaml $rundir
ln -s `pwd`/datafiles $rundir/datafiles

cd $rundir
${xfitter} >& xfitter.txt
cd - > /dev/null

flagBAD=0
# check chi2 in Results.txt ("After minimisation ...")
if [ $COPYRESULTS -eq 0 ]; then
  # some tests do not call 'fcn 3' and there is no chi2 stored in Results.txt
  grep  'After' ${EXAMPLEDIR}/Results.txt > temp/def.txt
  exitcode=$?
  if [ $exitcode = 0 ]; then
    grep  'After' $rundir/output/Results.txt > temp/out.txt
    exitcode=$?
    if [ $exitcode = 0 ]; then
      cat temp/out.txt
      diff temp/out.txt temp/def.txt 
      exitcode=$?
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
    else
        echo "========================================"
        echo -e "$FAILED no chi2 calculated: check temp/$TESTNAME/xfitter.txt"
        echo "========================================"
    fi
  fi
  rm -f temp/out.txt temp/def.txt
fi

# check all output files
if [ $COPYRESULTS -eq 0 ]; then
  echo "Checking all output files ..."
  for targetfile in `find ${EXAMPLEDIR} -type f`; do
    file=`echo ${targetfile} | sed -e s@${EXAMPLEDIR}@${rundir}/output@`
    if [ ! -f $file ]; then
      echo "No expected file $file"
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
  if [ $flagBAD == 0 ]; then
    echo -e "Everything is $OK"
  else
    echo -e "Something $FAILED: see above for details"
  fi
  echo "========================================"
else
  rm -rf $EXAMPLEDIR
  cp -r $rundir/output $EXAMPLEDIR
  echo "Output copied"
  echo "========================================"
fi

exit $flagBAD
