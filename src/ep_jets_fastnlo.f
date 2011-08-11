      Subroutine GetJetsFastNLOXsection(IDataSet)
C---------------------------------------------------------------------------
C
C  Created 01/08/2011.  Calculate ep jets cross sections
C
C---------------------------------------------------------------------------
      implicit none
      include 'ntot.inc'
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'

      integer IDataSet
      integer NFmax
      parameter(NFmax=100)
      integer i, idx
      double precision XSec(NFmax)

C  Get cross sections
      
      call fastnlocalc(IDataSet, XSec)
c      stop

      do i=1,NDATAPOINTS(IDataSet)
         
         idx =  DATASETIDX(IDataSet,i)
                  
         THEO(idx) =  XSec(i)

c         print *,DATEN(idx),THEO(idx)
      enddo

C      stop

      end
