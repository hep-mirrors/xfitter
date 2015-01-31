      Subroutine GetDYCCXsection(IDataSet)
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
      integer IDataSet
C-------------------------------------------------
      if (DATASETTheoryType(IDataSet).eq.'kfactor') then
         call GetDYCCXsection_kfactor(IDataSet)
      elseif (DATASETTheoryType(IDataSet).eq.'applgrid') then
         call GetDYXsection_applgrid(IDataSet)         
      endif
      end


      Subroutine GetDYNCXsection(IDataSet)
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
      integer IDataSet
C-------------------------------------------------
      if (DATASETTheoryType(IDataSet).eq.'kfactor') then
         call GetDYNCXsection_kfactor(IDataSet)
      elseif (DATASETTheoryType(IDataSet).eq.'applgrid') then
         call GetDYXsection_applgrid(IDataSet)         
      endif
      end

      
      

      Subroutine GetDYXsection_applgrid(IDataSet)         
C---------------------------------------------------
C
C Calculate DY CC cross section using APPLGRID interface.
C
C---------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "scales.inc"

      integer IDataSet      ! data set index
      integer NPMax         ! max. number of DY points
      parameter (NPMax = 500)
      double precision XSec(NPMax) ! applgrid convolution result
      double precision XSecPlus(NPMax) ! applgrid convolution result                                              
      double precision XSecMinus(NPMax) ! applgrid convolution result                                             

      integer i,idx,idxUnit
      double precision TheoryUnit  ! scale factor for theory to bring to data units.

      integer GetInfoIndex         ! function thet returns index of an information string.

      integer nkfact, idxKFact

      integer GetKFactIndex

      double precision evalTotXS
      integer idxBinEta1,idxBinEta2
      double precision BinSize
      logical LAsymmetry
      integer GetBinIndex
      
C>18/01/2013----------------------
      integer idxNorm
      double precision CSpred_tot, aScaleNorm

C----------------------------------------------------      
      if (NDATAPOINTS(IDataSet).gt.NPmax) then
         print *,'ERROR IN GetDYXsection_applgrid'
         print *,'INCREASE NPmax to ',NDATAPOINTS(IDataSet)
         call HF_stop
      endif

C Check if normalisation is needed
      idxNorm = GetInfoIndex(IDataSet,'Normalised')
      if (idxNorm.gt.0) then
         aScaleNorm = DATASETInfo(idxNorm,IDataSet)
      endif

C     Check type of the data                                                                                      

      idx = GetInfoIndex(IDataSet,'asymmetry')
      if (idx.gt.0) then
         LAsymmetry = 
     $     DATASETInfo(idx,IDataSet).gt.0
      else
         LAsymmetry = .false.
      endif

      if (LAsymmetry) then
C     Asymmetry case                                                                                              
         call ag_convolute( DATASETApplgridTheoryIndex(1,IDataSet),
     $        DataSetIOrder(IDataSet),
     $        DataSetMuR(IDataSet),
     $        DataSetMuF(IDataSet),
     $        XSecPlus)
         
         call ag_convolute( DATASETApplgridTheoryIndex(2,IDataSet),
     $        DataSetIOrder(IDataSet),
     $        DataSetMuR(IDataSet),
     $        DataSetMuF(IDataSet),
     $        XSecMinus)
         
C     check if we have to divide APPLGRID prediction to convert units to data units:                              
         idxUnit = GetInfoIndex(IDataSet,'theoryunit')
         if (idxUnit.gt.0) then
            Theoryunit = DATASETInfo(idxUnit,IDataSet)
         else
            Theoryunit = 1.
         endif
      else
         call ag_convolute( DATASETTheoryIndex(IDataSet),
     $        DataSetIOrder(IDataSet),
     $        DataSetMuR(IDataSet),
     $        DataSetMuF(IDataSet),
     $        XSec)

C     check if we have to divide APPLGRID prediction to convert units to data units:
         idxUnit = GetInfoIndex(IDataSet,'theoryunit')
         if (idxUnit.gt.0) then
            Theoryunit = DATASETInfo(idxUnit,IDataSet)
         else
            Theoryunit = 1.
         endif
      endif
      

C Normalize if the field 'Normalised' takes place in the dataset
      if (idxNorm.gt.0.and.aScaleNorm.gt.0) then
            CSpred_tot = 0
C     Calculate the bin size:
         idxBinEta1 = GetBinIndex(IDataSet,'y1')
         idxBinEta2 = GetBinIndex(IDataSet,'y2')
C     Calculate the total Xsec:
         do i=1,NDATAPOINTS(IDataSet)
            idx =  DATASETIDX(IDataSet,i)
            BinSize = AbstractBins(idxBinEta2,idx)-AbstractBins(idxBinEta1,idx)
            CSpred_tot = CSpred_tot + (XSec(i)*BinSize / TheoryUnit)
         enddo
       endif


C     perhaps we need to multiply prediction by EW or NNLO or both k-factors

      nkfact = DATASETNKfactors(IDataSet)

      do i=1, NDATAPOINTS(IDataSet)
C     Bin Size                                                                                                    
         idx =  DATASETIDX(IDataSet,i)

         if (LAsymmetry) then
           THEO(idx) =
     $           (XSecPlus(i)-XSecMinus(i))/(XSecPlus(i)+XSecMinus(i))
C          print*,'prediction for theory: ', THEO(idx),'--- index: ',i                                            
         else if(idxNorm.gt.0.and.aScaleNorm.gt.0) then
           THEO(idx) = XSec(i) / (TheoryUnit *CSpred_tot) /aScaleNorm
         else
           THEO(idx) = XSec(i) / TheoryUnit     
         endif

         !> Check if k-factor is required, multiply prediction by all of them.                                    
         do idxkfact=1,nkfact
            Theo(idx) = Theo(idx)*kfactors(idxkfact,idx)
         enddo
cc         print *,'hady', idx, THEO(idx), DATEN(idx),TheoryUnit                                                  
      enddo



      end

      Subroutine GetDYCCXsection_kfactor(IDataSet)
C-----------------------------------------------------------
C
C Created by SG, 26/05/2011
C Calculate DY W+, W- and asymmetry cross sections
C
C------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
      integer IDataSet
C-------------------------------------------
      integer NPmax
      parameter(NPmax=100)
      

      double precision wmp_bsigs(NPmax),BinSize


      integer i,idx
      integer idxKfactWplus, idxKfactWminus
      logical LAsymmetry,LWplus
      integer idxBinEta1,idxBinEta2
C Functions:
      integer GetInfoIndex
      integer GetKFactIndex
      integer GetBinIndex
C------------------------------------------------------------

      if (NDATAPOINTS(IDataSet)*2.gt.NPmax) then
         print *,'ERROR IN GetDYCCXsection'
         print *,'INCREASE NPmax to ',NDATAPOINTS(IDataSet)*2
         call HF_stop
      endif


c      if (IFlagFCN.eq.1) then
c         call initDYtmp(IDataSet)
c      endif


      call dy_get_res(IDataSet, wmp_bsigs)

      idxKfactWplus  = GetKFactIndex(IDataSet,'Wplus')
      idxKfactWminus = GetKFactIndex(IDataSet,'Wminus') 

C Apply k-factors:
      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
         if (idxKfactWminus.gt.0) then
            wmp_bsigs(2*i-1) = wmp_bsigs(2*i-1)*Kfactors(idxKfactWminus,idx)
         endif
         if (idxKfactWplus.gt.0) then
            wmp_bsigs(2*i) = wmp_bsigs(2*i)*Kfactors(idxKfactWplus ,idx)
         endif
      enddo

C Check type of the data
      LAsymmetry =  
     $     DATASETInfo( GetInfoIndex(IDataSet,'asymmetry'), IDataSet).gt.0

      if (.not. LAsymmetry) then
         LWplus = DATASETInfo( GetInfoIndex(IDataSet,'e charge'), IDataSet).gt.0
C Need also bin sizes:
         idxBinEta1 = GetBinIndex(IDataSet,'eta1')
         idxBinEta2 = GetBinIndex(IDataSet,'eta2')
      endif

      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
         if (LAsymmetry) then
            THEO(idx) = (wmp_bsigs(2*i)-wmp_bsigs(2*i-1) )/(wmp_bsigs(2*i)+wmp_bsigs(2*i-1))
         else
            BinSize = AbstractBins(idxBinEta2,idx)-AbstractBins(idxBinEta1,idx)
            if (LWplus) then
               THEO(idx) = wmp_bsigs(2*i)/BinSize
            else
               THEO(idx) = wmp_bsigs(2*i-1)/BinSize         
            endif
         endif
      enddo
      end


      subroutine  GetDYNCXsection_kfactor(IDataSet)
C-----------------------------------------------------------
C
C Created by SG, 26/05/2011
C C&P from CC by AS 06/07/2011
C The normalization of the predictions is added 15/11/2011
C
C------------------------------------------------------------
      implicit none
#include "ntot.inc"
#include "steering.inc"
#include "for_debug.inc"
#include "datasets.inc"
#include "indata.inc"
#include "theo.inc"
#include "fcn.inc"
      integer IDataSet
C-------------------------------------------
      integer NPmax
      parameter(NPmax=100)

      double precision z_bsigs(NPmax),BinSize

      integer i,idx
      integer idxKfactZ
      logical LAsymmetry,LZ
      integer idxBinY1,idxBinY2

C>15/11/2011----------------------
      integer idxNorm
      double precision CSpred_tot, aScaleNorm
C<15/11/2011----------------------

C Functions:
      integer GetInfoIndex
      integer GetKFactIndex
      integer GetBinIndex
C------------------------------------------------------------

      if (NDATAPOINTS(IDataSet)*2.gt.NPmax) then
         print *,'ERROR IN GetDYNCXsection'
         print *,'INCREASE NPmax to ',NDATAPOINTS(IDataSet)*2
         call HF_stop
      endif

C>15/11/2011----------------------
C Check if normalisation is needed
      idxNorm = GetInfoIndex(IDataSet,'Normalised')
      if (idxNorm.gt.0) then
         aScaleNorm = DATASETInfo(idxNorm,IDataSet)
      endif
C<15/11/2011----------------------

c      if (IFlagFCN.eq.1) then
c         call initDYtmp(IDataSet)
c      endif


      call dy_get_res(IDataSet, z_bsigs)
      idxKfactZ  = GetKFactIndex(IDataSet,'Z0')


C Apply k-factors:
      if (idxKfactZ.gt.0) then
         do i=1,NDATAPOINTS(IDataSet)
            idx =  DATASETIDX(IDataSet,i)
            z_bsigs(i) = z_bsigs(i)*Kfactors(idxKfactZ,idx)
         enddo
      else
         print '(''GetDYNCXsection_kfactor: no kfactor found'')'
      endif


C Need also bin sizes:
      idxBinY1 = GetBinIndex(IDataSet,'y1')
      idxBinY2 = GetBinIndex(IDataSet,'y2')

C>15/11/2011----------------------
C Normalize if the field 'Normalised' takes place in the dataset
      if (idxNorm.gt.0.and.aScaleNorm.gt.0) then
         CSpred_tot = 0
C Calculate sum:
         do i=1,NDATAPOINTS(IDataSet)
            CSpred_tot = CSpred_tot + z_bsigs(i)
         enddo
C Normalise:
         do i=1,NDATAPOINTS(IDataSet)
            z_bsigs(i) = z_bsigs(i)/CSpred_tot/aScaleNorm
         enddo
      endif
C<15/11/2011----------------------


      do i=1,NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
c         if (LAsymmetry) then
c            THEO(idx) = (wmp_bsigs(2*i)-wmp_bsigs(2*i-1) )/(wmp_bsigs(2*i)+wmp_bsigs(2*i-1))
c         else
            BinSize = AbstractBins(idxBinY2,idx)-AbstractBins(idxBinY1,idx)
            THEO(idx) = z_bsigs(i)/BinSize
c         endif
      enddo
      end


