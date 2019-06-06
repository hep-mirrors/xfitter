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

rm -rf temp
mkdir temp

echo "========================================"
echo "Running checks"
echo "========================================"
echo "validation test: "
echo "PASS if code runs properly DIS code and produces HERAPDF2.0"
echo "FAIL if code fails to reproduce HERAPDF2.0"
echo "========================================"

INPUTDIR='input_steering_22'
EXAMPLEDIR='examples_22'

cp ${INPUTDIR}/steering.txt.def steering.txt
cp ${INPUTDIR}/parameters.yaml.def parameters.yaml
cp ${INPUTDIR}/constants.yaml.def constants.yaml

bin/xfitter >/dev/null

grep  'After' output/Results.txt > temp/out.txt
grep  'After' ${EXAMPLEDIR}/output/Results.txt > temp/def.txt

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
  out=`checkFile $file ${EXAMPLEDIR}/${file}`
  echo $out
  echo $out | grep FAILED > /dev/null
  if [ $exitcode = 1 ]; then
    flagAllFine=0
  fi
done
echo "========================================"
if [ $flagAllFine = 1 ]; then
  echo "Everything is fine"
else
  echo "Something failed: see above for details"
fi
echo "========================================"

rm -rf temp

if [ $flagAllFine = 0 ]; then
  exit 1
fi

