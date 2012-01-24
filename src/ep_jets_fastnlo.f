      Subroutine GetJetsFastNLOXsection(IDataSet, UseNormalisation)
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
      logical UseNormalisation
      integer NFmax
      parameter(NFmax=100)
      integer i, idx, idxNPCorr, idxZ0Corr
      double precision XSec(NFmax), XNP, XZ0

C Functions:
      integer GetBinIndex

C  Get cross sections
      
      call fastnlocalc(IDataSet, XSec)
      idxNPCorr = GetBinIndex(IDataSet,'NPCorr')
      idxZ0Corr = GetBinIndex(IDataSet,'Z0Corr')

      do i=1,NDATAPOINTS(IDataSet)
         
         idx =  DATASETIDX(IDataSet,i)
         
         if (idxNPCorr.gt.0) then
            XNP = AbstractBins(idxNPCorr,idx)
         else
            XNP = 1.0
         endif

         if (idxZ0Corr.gt.0) then
            XZ0 = AbstractBins(idxZ0Corr,idx)
         else
            XZ0 = 1.0
         endif

         if(UseNormalisation) then
            THEO(idx) =  XSec(i) * XNP * XZ0 / THEO(idx)
         else
            THEO(idx) =  XSec(i) * XNP * XZ0
         endif
c         print *,DATEN(idx),THEO(idx), XNP, XZ0
      enddo

      end


      Subroutine GetQcdnumPdfset(pdfset)
C---------------------------------------------------------------------------
C
C  Created 05/12/2011.  Get pdfset to be used by qcdnum
C
C---------------------------------------------------------------------------
      implicit none
      include 'steering.inc'
      integer pdfset
      pdfset = IPDFSET
      end
