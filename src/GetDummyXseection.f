      Subroutine GetDummyXsection(IDataSet)



      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"

      integer i,idx
      integer IDataSet      ! data set index

      do i=1, NDATAPOINTS(IDataSet)

         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) = DATEN(idx)
      enddo
      end
