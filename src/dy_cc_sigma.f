      Subroutine GetDYCCXsection(IDataSet)
      implicit none
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
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
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
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
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'

      integer IDataSet      ! data set index
      integer NPMax         ! max. number of DY points
      parameter (NPMax = 200)
      double precision XSec(NPMax) ! applgrid convolution result

      integer i,idx,idxUnit
      double precision TheoryUnit  ! scale factor for theory to bring to data units.

      integer GetInfoIndex         ! function thet returns index of an information string.
C----------------------------------------------------      
      call ag_convolute( DATASETTheoryIndex(IDataSet),XSec)

C check if we have to divide APPLGRID prediction to convert units to data units:
      idxUnit = GetInfoIndex(IDataSet,'theoryunit')
      if (idxUnit.gt.0) then
         Theoryunit = DATASETInfo(idxUnit,IDataSet)
      else
         Theoryunit = 1.
      endif

      do i=1, NDATAPOINTS(IDataSet)
         idx =  DATASETIDX(IDataSet,i)
         THEO(idx) = XSec(i) / TheoryUnit

c         print *,'hady', idx, THEO(idx), DATEN(idx)
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
      include 'steering.inc'
      include 'for_debug.inc'
      include 'datasets.inc'
      include 'ntot.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
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
         stop
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
      implicit none
      integer IDataSet
C--------------------------------------------------
      print *,'Call dummy GetDYNCXsection_kfactor'
      end


c      subroutine Initdytmp(IDataSet)
cC------------------------------------------------------------
cC
cC Initialise tables for DY process
cC
cC------------------------------------------------------------
c      implicit none
c      include 'steering.inc'
c      include 'for_debug.inc'
c      include 'datasets.inc'
c      include 'ntot.inc'
c      include 'indata.inc'
c 
c
c      integer IDataSet
c      integer GetBinIndex                                                                                                                                    
c      double precision dy_mass(2)                                                                                                                            
c      double precision dy_y(2)
c      double precision pt_cut
c      integer NPmax
c      parameter(NPmax=100)
c      double precision eb(Npmax+1)
c
c      integer idx, idxEta1, idxEta2,i 
c
cC----------------------------------------------------------
c
c      if (NDATAPOINTS(IDataSet).gt.NPmax) then
c         print *,'ERROR IN InitDYCCXsectionDataset'
c         print *,'INCREASE NPMax to ',NDATAPOINTS(IDataSet)                                                                                                  
c         stop
c      endif
c
c
c      pt_cut=25.
c      dy_mass(1) = 1.
c      dy_mass(2) = 7000.
c      
c      dy_y(1) = -10.
c      dy_y(2) =  10.
c
cC Get indicies:
c      idxEta1 = GetBinIndex(IDataSet,'eta1')
c      idxEta2 = GetBinIndex(IDataSet,'eta2')
c
c      if (idxEta1.eq.0 .or. idxEta2.eq.0) then
c         print 
c     $        '(''ERROR in GetDYCCXsection, can not find bin index for Eta1, Eta2'',2i6)'
c     $        ,idxEta1,idxEta2
c         stop
c      endif
c
cC Define bins:
c      do i=1,NDATAPOINTS(IDataSet) 
c         idx =  DATASETIDX(IDataSet,i)
c         eb(i)   =  AbstractBins(idxEta1,idx)
c         eb(i+1) =  AbstractBins(idxEta2,idx)
c      enddo
c
c      print *,'Initialise DY calculations for dataset', IDataSet
c      call dy_create_calc(IDataSet, 1, 7000d0, 'W'//char(0), dy_mass, dy_y, pt_cut)
c      call dy_set_bins(IDataSet,'eta'//char(0), NDATAPOINTS(IDataSet), eb)
c
c      end

