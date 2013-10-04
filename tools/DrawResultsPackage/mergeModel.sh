#!/bin/bash

PARAMNAMES="fs_ mb_ mc_ q2m_"
PARAMVARS="m p"
Q2NBINS=6
DIRBASE="../base"

PARNUM=1

if [ -d "summ" ];
 then
  rm -r summ
fi
mkdir summ

cp $DIRBASE/parsout_0 $DIRBASE/fittedresults.txt $DIRBASE/lhapdf.block.txt $DIRBASE/pulls.first.txt $DIRBASE/pulls.last.txt $DIRBASE/Results.txt $DIRBASE/pdfs_q2val_*.txt summ/

for PARN in $PARAMNAMES; do
  if [ "$PARNUM" -lt "10" ];
   then
    PARNUMST="0${PARNUM}"
   else
    PARNUMST="${PARNUM}"
  fi
  for PARV in $PARAMVARS; do
    for (( Q2B=1; Q2B<=Q2NBINS; Q2B++ )); do
      if [ "$Q2B" -lt "10" ];
       then
        Q2BST="0${Q2B}"
       else
        Q2BST="${Q2B}"
      fi
      cp ${PARN}${PARV}/pdfs_q2val_${Q2BST}.txt summ/pdfs_q2val_m${PARNUMST}${PARV}_${Q2BST}.txt
    done
  done
  let "PARNUM += 1"
done
