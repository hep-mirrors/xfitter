      Subroutine GetHathorXsection(IDataSet)
C---------------------------------------------------------------------------
C
C  Created 09/11/2011.  Calculate ttbar cross sections
C
C---------------------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"

      integer IDataSet
      double precision XSec

      call hathorcalc(IDataSet, XSec)

      THEO(DATASETIDX(IDataSet,1)) = Xsec

      end
