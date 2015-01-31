      Subroutine GetJetsPPApplGrid(IDataSet)
C---------------------------------------------------------------------------
C
C  Created 03/06/2011.  Convolution of apple grid for PP jet cross sections
C
C---------------------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "scales.inc"

      integer IDataSet
      integer NPMax,NEtaMax
      parameter (NPMax=200)
      parameter (NEtaMax=20)

C Convolution results:
      double precision XSec(NPMax,NEtaMax),XS,XNP,BS

      integer i,idx, idxNPCorr, idxYBinSize
      integer idxEta,iEtaMin, iEtaMax,iEta

C Functions:
      integer GetBinIndex
      integer GetInfoIndex


C-------------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetJetsPPApplGrid'
         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif

      idxEta = GetBinIndex(IDataSet,'EtaBinNumber')

C
C Arrange loop over eta bins:
C
      if (idxEta.gt.0) then
         iEtaMax = 0
         iEtaMin = 1000
         do i=1,NDATAPOINTS(IDataSet)
            idx =  DATASETIDX(IDataSet,i)
            iEta = int(AbstractBins(idxEta,idx)+0.1)
            if (iEta.gt.iEtaMax) iEtaMax = iEta
            if (iEta.lt.iEtaMin) iEtaMin = iEta
         enddo
      else
         iEtaMin = 1
         iEtaMax = 1
      endif
      
      if (iEtaMax.gt.NEtaMax) then
         print *,'Error in GetJetsPPApplGrid'
         print *,'INCREASE NEtaMax'
         call HF_stop
      endif


      do iEta=iEtaMin,iEtaMax
C Convolution:
         call ag_convolute( DATASETTheoryIndex(IDataSet)+iEta-1,
     $        DataSetIOrder(IDataSet),
     $        DataSetMuR(IDataSet),
     $        DataSetMuF(IDataSet),
     $        XSec(1,iEta))
      enddo

C Check NP correction:
      idxNPCorr = GetBinIndex(IDataSet,'NPCorr')
C Check Bin size in rapidity:
      idxYBinSize = GetBinIndex(IDataSet,'YBinSize')

      do i=1,NDATAPOINTS(IDataSet)

         idx =  DATASETIDX(IDataSet,i)
         
         if (idxEta.gt.0) then         
            iEta = int(AbstractBins(idxEta,idx)+0.1)
         else
            iEta = 1
         endif
         
         XS  =  XSec( IndexTheoryBin(idx),iEta )
         
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

C      call HF_stop

      end
