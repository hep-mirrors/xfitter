#!/bin/sh
#
#   Turn on Verbose Mode
#set -x

rm -rf temp
mkdir temp

echo "========================================"
echo "Running checks"
echo "========================================"
echo "validation test: "
echo "PASS if code runs properly DIS code and produces HERAPDF2.0"
echo "FAIL if code fails to reproduce HERAPDF2.0"
echo "========================================"


cp input_steering/minuit.in.txt.def minuit.in.txt
cp input_steering/steering.txt.def steering.txt

bin/xfitter >/dev/null

grep  'After' output/Results.txt > temp/out.txt
grep  'After' examples/output/Results.txt > temp/def.txt

diff temp/out.txt temp/def.txt 

exitcode=$?

if [ $exitcode = 0 ]; then
echo "========================================"
echo "Check is fine"
echo "========================================"
else
echo "========================================"
echo -e "Failed validation with default steering"
echo "========================================"
exit 1
fi 
echo "========================================"

rm -rf temp


