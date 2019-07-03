#!/bin/bash
#
# Collection of all tests (slow tests might be commented out)
#
# If output for some tests needs to be updated, add '--copy', e.g.:
# ./tools/check22-all.sh --copy
#

COPY=''
if [ "$1" = "--copy" ]; then
  COPY=' --copy'
fi

echo "====================================================="
echo "Running all tests ... output will be stored in temp/"
echo "====================================================="

rm -rf temp/
flagAllFine=1

# chi2 iteration NNLO RTOPT QCDNUM (default) [HERAPDF2.0 arXiv:1506.06042 ]
# it tests also storing NNLO PDF from QCDNUM evolution in LHAPDF6 format
./tools/check22.sh defaultNNLO $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration NLO RTOPT QCDNUM [HERAPDF2.0 arXiv:1506.06042 ]
# it tests also storing NLO PDF from QCDNUM evolution in LHAPDF6 format
./tools/check22.sh defaultNLO $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration NLO FFABM QCDNUM [HERAPDF-HVQMASS arXiv:1804.01019]
./tools/check22.sh FFABM $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration FONLL APFEL
# it tests also storing PDF from APFEL evolution in LHAPDF6 format
./tools/check22.sh FONLL $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration for ALL data using NLO RTOPT LHAPDF=CT10nlo
# it tests also storing NNLO PDF from LHAPDF evolution in LHAPDF6 format
./tools/check22.sh ALLDATA $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration for AFB pseudodata, checks AFB (LO) and APPLgrid (NLO) reactions
./tools/check22.sh AFB $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration for LHeC pseudodata
./tools/check22.sh LHeC $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration for PROSA analysis 1503.04581 (absolute variant)
./tools/check22.sh HVQMNR-abs $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration for PROSA analysis 1503.04581 (normalised variant)
./tools/check22.sh HVQMNR-norm $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# chi2 iteration for small x resummation xfitter analysis 1802.00064 
./tools/check22.sh SmallxResummation $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# profiler alpha_s
./tools/check22.sh profilerAs $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# profiler PDF eigenvectors (LHAPDF)
./tools/check22.sh profilerLHAPDF $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

# full fit ZMVFNS NNLO QCDNUM, with error bands (takes ~ 10 min)
# it tests also storing PDF eigenvectors in LHAPDF6 format
./tools/check22.sh ZMVFNS-fit $COPY
if [ `echo $?` -ne 0 ]; then flagAllFine=0; fi

echo "====================================================="
if [ $flagAllFine -eq 1 ]; then
  if [ -z $COPY ]; then
    echo " -> All tests are fine"
  else
    echo " -> All output copied"
  fi
else
  echo " -> Something failed: see above for details (logs are in /temp)"
fi
echo "====================================================="

if [ $flagAllFine = 0 ]; then
  exit 1
fi
