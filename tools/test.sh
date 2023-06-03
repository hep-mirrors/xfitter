#!/bin/bash

# list of tests to omit (if commented out, no tests are omitted)
#omitTests=('ZMVFNS-fit' 'profilerLHAPDF') # these are two slow tests, skipping them will save ~15min
#omitTests=('ceresZMVFNSfastChi2' 'chi2scanMTOP')
omitTests=('profilerCIJET' 'ZPT' 'ZMVFNS-fit' 'scanmin' 'CERES-fit' 'CERES-parallel' 'CERES-Chebyschev') 

install_dir=$(pwd)
# xfitter binary
xfitter=$install_dir/'bin/xfitter'

# xfitter libs
export LD_LIBRARY_PATH=$install_dir/lib:$install_dir/lib/xfitter:$LD_LIBRARY_PATH

# log file for xfitter output
xflogfile='xfitter.log'

# log file for test output
testlogfile='test.log'

# output messages
FAILED="\e[31m\e[1mFAIL\e[0m" # bold red
PASSED="\e[32m\e[1mPASS\e[0m" # bold green

# function to check if element $1 is in array $2
# code from https://stackoverflow.com/questions/3685970/check-if-a-bash-array-contains-a-value
containsElement () {
  local e match="$1"
  shift
  for e; do [[ "$e" == "$match" ]] && return 0; done
  return 1
}

# tolerate small differences in terms of maximum fraction of different lines
# some PDF values in files output/pdfs_q2val_ and LHAPDF output files could differ,
# but typically there are very few different lines output/pdfs_q2val_* and */*_*.dat (LHAPDF files)
# differences are allowed only in files 
function tolerateDiff()
{
  if [[ $1 != *"output/pdfs_q2val_"* ]] && [[ $1 != *".dat" ]]; then
    return 1;
  fi
  maxFractionDiff=$3
  # make sure two files have the same number of lines
  nl1=`cat $1 | wc -l`
  nl2=`cat $2 | wc -l`
  if [ $nl1 -ne $nl2 ]; then 
    return 1;
  fi
  # calculate which fraction of lines is different
  nldiff=`$diff $1 $2 | grep "^>" | wc -l`
  fractionDiff=`echo $nldiff/$nl1 | bc -l`
  #echo "tolerateDiff $1 $fractionDiff"
  status=`echo "$fractionDiff>$maxFractionDiff" | bc -l` # will be 1 if difference is too large and test not passed
  return $status
}

# function to check if $1 and $2 are identical
checkFile()
{
  printf "$diff $1 $2 ... "
  if [[ $diff == *'numdiff'* ]]; then
      $diff -s ' \t\n,' $1 $2 > /dev/null
  else
      $diff  $1 $2 > /dev/null
  fi
  exitcode=$?
  
  #if [ $exitcode = 1 ]; then
  #  # check if we can tolerate small differences in some tests
  #  # allow not more than 2% of different lines in some output files of ZMVFNS-fit test
  #  # (see tolerateDiff for details which files could be different)
  #  if [[ $1 = *"ZMVFNS-fit"* ]]; then
  #    tolerateDiff $1 $2 0.02
  #    exitcode=$?
  #  fi
  #fi
  
  if [ $exitcode = 0 ]; then
    echo "PASSED"
  else
    echo "FAILED"
  fi 
}

runTest()
{
  # test status
  flagBAD=0

  # TESTNAME is the name of directory in examples/
  TESTNAME=$1
  TESTNAMEDEF='defaultNNLO'
  if [ -z $TESTNAME ]; then
    TESTNAME=$TESTNAMEDEF
  fi
  # directory where to run
  rundir=$2
  # Optionally copy results and make them reference
  COPYRESULTS=$3
  #rm -rf $rundir
  #mkdir -p $rundir
  
  # diff command: by defult use 'diff', but if available use 'numdiff' which tolerate small differences
  diff='diff'
  # we want this numdiff: https://www.nongnu.org/numdiff/
  # check "numdiff -v" command: avoid using some other "numdiff" installed on some systems (like presently on naf-xfitter.desy.de) which does not recognize "-v"
  if [ `numdiff -v >& /dev/null; echo $?` == "0" ]; then
    # we have numdiff and we wil use it with tolerance 1e-4 for either absolute or relative differneces between numbers
    diff='numdiff -a 1e-3 -r 1e-3'
  fi

  echo "========================================"
  echo "Running check: $TESTNAME"
  echo "========================================"
  if [ $COPYRESULTS -eq 0 ]; then
    echo "This is validation test:"
    echo "OK if code runs properly"
    echo "FAILED if code fails to reproduce expected results"
    echo "========================================"
  fi

  INPUTDIR="examples/$TESTNAME"
  EXAMPLEDIR="examples/$TESTNAME/output"
  if [ $COPYRESULTS -eq 1 ]; then
    echo "Results will be stored as reference in $EXAMPLEDIR"
  else
    if [ ! -d $EXAMPLEDIR ]; then
      echo "Warning: no reference output directory -> test will be considered FAILED"
      flagBAD=1
    fi
  fi
  echo "Running in temp/$TESTNAME"
  echo "Log file stored in temp/$TESTNAME/$xflogfile"

  if [ ! -d $INPUTDIR ]; then
    echo "Failed to find input files for test \"$TESTNAME\""
    echo "expected directory $INPUTDIR"
    return 1
  fi
  echo "Using input files from ${INPUTDIR}"
  echo "========================================"

  cp ${INPUTDIR}/steering.txt $rundir
  cp ${INPUTDIR}/parameters.yaml $rundir
  cp ${INPUTDIR}/constants.yaml $rundir
  # also copy any .dat files
  cp ${INPUTDIR}/*.dat $rundir
  ln -s `pwd`/datafiles $rundir/datafiles

  cd $rundir
  ${xfitter} >& ${xflogfile}
  cd - > /dev/null

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
        $diff temp/out.txt temp/def.txt 
        exitcode=$?
        if [ $exitcode = 0 ]; then
          echo "========================================"
          echo "Check of chi^2 is PASSED"
          echo "========================================"
        else
          echo "========================================"
          echo "FAILED validation with default steering"
          echo "========================================"
          flagBAD=1
        fi
      else
          echo "========================================"
          echo "FAILED no chi2 calculated: check temp/$TESTNAME/$xflogfile"
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
      echo "Everything is PASSED"
    else
      echo "Something FAILED: see above for details"
    fi
    echo "========================================"
  else
    rm -rf $EXAMPLEDIR
    cp -r $rundir/output $EXAMPLEDIR
    echo "Output copied"
    echo "========================================"
  fi

  return $flagBAD
}

##########################
# the script starts here #
##########################
if [ $# -ne 0 ] && [ "${@: -1}" = "--help" ]; then
  echo "Usage: test.sh <TEST1> <TEST2> ... [OPTION]"
  echo "OPTION could be:"
  echo "  --copy to copy results and make them reference"
  echo "  --help to see this message"
  echo "If no test names are provided, all tests in directory examples/ will run, except those specified in 'omitTests' list"
  exit 1
fi

COPY=0
if [ $# -ne 0 ] && [ "${@: -1}" = "--copy" ]; then
  COPY=1
  echo "==========================================================================="
  echo "Running in COPY mode: output of tests will be copied and saved as reference"
  echo "==========================================================================="
fi

listOfTests=""
for arg in "$@"; do
  if [ "${arg:0:2}" != "--" ]; then
    listOfTests="$listOfTests $arg"
  fi
done
if [ -z "$listOfTests" ]; then
  for dir in `ls -1d examples/*`; do
    # skip possible files: only directories contain tests
    if [ ! -d $dir ]; then continue; fi
    dir=`basename $dir`;
    containsElement $dir "${omitTests[@]}"
    if [ `echo $?` -eq 1 ]; then
      listOfTests="$listOfTests $dir"
    fi
  done
fi
if [ -z "$listOfTests" ]; then
  echo "No tests to run"
fi

testsPassed=0
testsFailed=0
for arg in `echo $listOfTests`; do
  dir=temp/$arg
  rm -rf $dir
  mkdir -p $dir
  log=$dir/$testlogfile
  printf "Testing $arg ... "
  rm -rf 
  runTest $arg $dir $COPY >& $log
  exitcode=$?
  if [ $COPY -eq 1 ]; then
    if [ $exitcode == 0 ]; then
      echo copied [details in $log]
      testsPassed=$[$testsPassed+1]
    else
      echo not copied [details in $log]
      testsFailed=$[$testsFailed+1]
    fi
  else
    if [ $exitcode == 0 ]; then
      echo -e $PASSED [details in $log]
      testsPassed=$[$testsPassed+1]
    else
      echo -e $FAILED [details in $log]
      testsFailed=$[$testsFailed+1]
    fi
  fi
done

if [ $COPY -eq 1 ]; then
  if [ $testsPassed -gt 0 ]; then
    echo -e "-> $testsPassed test(s) copied"
  fi
  if [ $testsFailed -gt 0 ]; then
    echo -e "-> $testsFailed test(s) not copied"
  fi
else
  if [ $testsPassed -gt 0 ]; then
    echo -e "-> $testsPassed test(s) $PASSED"
  fi
  if [ $testsFailed -gt 0 ]; then
    echo -e "-> $testsFailed test(s) $FAILED"
  fi
fi

exit $testsFailed
