      Subroutine GetJetsPPApplGrid(IDataSet)
C----------------------------------------------------------------
C
C  Created 03/06/2011.  Convolution of apple grid
C
C----------------------------------------------------------------
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'

      integer IDataSet
      integer NPMax
      parameter (NPMax=100)

C Convolution results:
      double precision XSec(NPMax),XS,XNP,BS

      integer i,idx, idxNPCorr, idxYBinSize

C Functions:
      integer GetBinIndex
      integer GetInfoIndex


C-------------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetJetsPPApplGrid'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         stop
      endif

C Convolution:
      call ag_convolute( DATASETTheoryIndex(IDataSet), XSec)

C Check NP correction:
      idxNPCorr = GetBinIndex(IDataSet,'NPCorr')
C Check Bin size in rapidity:
      idxYBinSize = GetBinIndex(IDataSet,'YBinSize')

      do i=1,NDATAPOINTS(IDataSet)

         idx =  DATASETIDX(IDataSet,i)
         
         XS  =  XSec( IndexTheoryBin(idx) )
         
         if (idxNPCorr.gt.0) then
C apply extra NP correction:
            XNP = AbstractBins(idxNPCorr,idx)
         else
            XNP = 1.0
         endif
         

         if (idxYBinSize.gt.0) then
C apply extra binsize correction:
            BS = AbstractBins(idxYBinSize,idx)
         else
            BS = 1.
         endif

c         print *,'check',i,IndexTheoryBin(idx),idxNPCorr,XS,XNP,BS

         THEO(idx) =  XS/BS*XNP

c         print *,DATEN(idx),THEO(idx)
         
      enddo

c      stop

      end
