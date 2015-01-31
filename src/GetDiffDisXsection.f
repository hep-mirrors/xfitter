cjt      Subroutine GetDiffDisXsection(IDataSet, XSecType)
	Subroutine GetDiffDisXsection(IDataSet)
C----------------------------------------------------------------
C
C  NC triple differential reduced cross section calculation 
C  for dataset IDataSet. Fills global array THEO.
C  Needs 'MX', 'Q2', 'xpom' columns in a data file and following CInfo fields:
C                'sqrt(S)'
C
C  \author Wojtek Slominski, Jagiellonian Univ., Physics Dept.
C  \date 2012-03
C  \copyright Creative Commons license CC-BY-NC 3.0 
C---------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
#include "couplings.inc"
#include "qcdnumhelper.inc"

cjt      character*(*) XSecType
      integer IDataSet
      integer idQ2, idMX, idxP, i,  idx
     
      
      double precision xPom(NPMaxDIS),MX(NPMaxDIS),Q2(NPMaxDIS)
      double precision XSec(NPMaxDIS)
cws      double precision Charge, polarity, alphaem_run, factor
cws      logical IsReduced
      double precision vars(3)

C Functions:
      integer GetBinIndex
      integer GetInfoIndex
      double precision DDISvalue

c H1qcdfunc
      integer ifirst
      data ifirst /1/
C---------------------------------------------------------


      if (NDATAPOINTS(IDataSet).gt.NPMaxDIS) then
         print *,'ERROR IN GetDiffDisXsection'
         print *,'INCREASE NPMaxDIS to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif

C
C Get indexes for Q2, xPom and MX bins:
C
      idQ2 = GetBinIndex(IDataSet,'Q2')
      idxP  = GetBinIndex(IDataSet,'xpom')
      idMX = GetBinIndex(IDataSet,'MX')


      if (idQ2.eq.0 .or. idxP.eq.0 .or. idMX.eq.0) then
         Return
      endif

C prepare bins:
      do i=1,NDATAPOINTS(IDataSet)
C
C Reference from the dataset to a global data index:
C
         idx =  DATASETIDX(IDataSet,i)
C
C Local xPom,MX,Q2 arrays, used for QCDNUM SF caclulations:
C
	 
         xPom(i)   = AbstractBins(idxP,idx)
         MX(i)     = AbstractBins(idMX,idx)
         Q2(i)     = AbstractBins(idQ2,idx)
      enddo


      do i=1,NDATAPOINTS(IDataSet)
        vars(1) = xPom(i)
        vars(2) = Q2(i)
        vars(3) = MX(i)
        XSec(i) = DDISvalue(IDataSet,vars)
      enddo


      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) =  XSec(i)
      enddo

      end
